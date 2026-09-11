import gzip
import logging
import os
import pickle
import tempfile
from heapq import heappush

try:
    from .pvalue_handler import PvalueHandler
except ImportError:  # Executed as the upstream standalone script.
    from pvalue_handler import PvalueHandler



valid_nucleotides = {"A", "T", "G", "C"}


class ReferenceSequenceFetchError(RuntimeError):
    """The reference sequence could not be retrieved from the configured FASTA."""


class ReferenceAlleleMismatchError(ValueError):
    """The supplied reference allele differs from a successfully fetched sequence."""


class InvalidAlleleError(ValueError):
    """An allele contains a character outside the supported DNA alphabet."""


def normalize_alleles(reference, start, stop, alleles):
    """Return a minimal, left-aligned representation within ``reference``.

    This is the two-allele subset of the normalization used by the upstream
    gwas2vcf adapter. Keeping it here avoids its abandoned compiled ``vgraph``
    dependency while retaining the same trim-then-left-shuffle behaviour.
    Coordinates are zero-based and half-open.
    """
    alleles = [str(allele).upper() for allele in alleles]
    if alleles[0] != reference[start:stop].upper():
        raise ValueError("Reference allele does not match the reference sequence")

    while all(alleles) and len({allele[-1] for allele in alleles}) == 1:
        alleles = [allele[:-1] for allele in alleles]
        stop -= 1
    while all(alleles) and len({allele[0] for allele in alleles}) == 1:
        alleles = [allele[1:] for allele in alleles]
        start += 1

    while any(not allele for allele in alleles) and start > 0:
        base = reference[start - 1].upper()
        alleles = [base + allele for allele in alleles]
        start -= 1
        stop -= 1
        if len({allele[-1] for allele in alleles}) != 1:
            break
        alleles = [allele[:-1] for allele in alleles]

    return start, stop, tuple(alleles)


class Gwas:
    def __init__(
        self,
        chrom,
        pos,
        ref,
        alt,
        b,
        se,
        nlog_pval,
        n,
        alt_freq,
        study_id,
        ncase,
        ncontrol,  # Added ncontrol
        neff,      # Added neff
        imp_info,
        imp_z,
        vcf_filter="PASS",
    ):
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.b = b
        self.se = se
        self.nlog_pval = nlog_pval
        self.n = n
        self.alt_freq = alt_freq
        self.study_id = study_id
        self.ncase = ncase
        self.ncontrol = ncontrol  # Store ncontrol
        self.neff = neff          # Store neff
        self.imp_info = imp_info
        self.imp_z = imp_z
        self.vcf_filter = vcf_filter

    def reverse_sign(self):
        ref_old = self.ref
        alt_old = self.alt
        self.ref = alt_old
        self.alt = ref_old
        self.b = self.b * -1
        if self.imp_z is not None:
            self.imp_z = self.imp_z * -1
        try:
            self.alt_freq = 1 - self.alt_freq
        except TypeError:
            self.alt_freq = None

    def check_reference_allele(self, fasta):
        start = self.pos - 1
        end = start + len(self.ref)
        try:
            contig_length = fasta.get_reference_length(self.chrom)
            if start < 0 or end > contig_length:
                raise ReferenceSequenceFetchError(
                    "Reference interval is outside the FASTA contig: "
                    f"{self.chrom}:{self.pos}-{self.pos + len(self.ref) - 1} "
                    f"(contig length {contig_length})"
                )
            fasta_ref_seq = fasta.fetch(
                reference=self.chrom,
                start=start,
                end=end,
            ).upper()
        except ReferenceSequenceFetchError:
            raise
        except Exception as exception_name:
            raise ReferenceSequenceFetchError(
                "Could not retrieve the reference sequence for "
                f"{self.chrom}:{self.pos}-{self.pos + len(self.ref) - 1}. "
                "Check that the chromosome naming, coordinates, FASTA, and FASTA "
                "index match the detected genome build."
            ) from exception_name

        if len(fasta_ref_seq) != len(self.ref):
            raise ReferenceSequenceFetchError(
                "The FASTA returned an incomplete reference sequence for "
                f"{self.chrom}:{self.pos}-{self.pos + len(self.ref) - 1}."
            )
        if self.ref != fasta_ref_seq:
            raise ReferenceAlleleMismatchError(
                f"Supplied REF {self.ref} does not match FASTA REF "
                f"{fasta_ref_seq} at {self.chrom}:{self.pos}"
            )

    def normalise(self, fasta, padding=100):
        if len(self.ref) < 2 and len(self.alt) < 2:
            return
        pos0 = self.pos - 1
        try:
            seq = fasta.fetch(
                reference=self.chrom, start=pos0 - padding, end=pos0 + padding
            ).upper()
        except Exception as exception_name:
            raise ReferenceSequenceFetchError(
                "Could not retrieve the normalization window for "
                f"{self.chrom}:{self.pos}. Check that the chromosome naming, "
                "coordinates, FASTA, and FASTA index match the detected genome build."
            ) from exception_name
        start, stop, alleles = normalize_alleles(
            seq, padding, padding + len(self.ref), (self.ref, self.alt)
        )
        self.ref = alleles[0]
        self.alt = alleles[1]
        self.pos = (pos0 - padding) + start + 1
        if len(self.ref) == 0 or len(self.alt) == 0:
            dist = (self.pos - 1) - pos0
            left_nucleotide = seq[(padding + dist) - 1 : (padding + dist)]
            self.ref = left_nucleotide + self.ref
            self.alt = left_nucleotide + self.alt
            self.pos = self.pos - 1

    def check_alleles_are_valid(self):
        if not self.alt or not self.ref:
            raise InvalidAlleleError("REF and ALT alleles must both be non-empty")
        for nucleotide in self.alt:
            if nucleotide not in valid_nucleotides:
                raise InvalidAlleleError(
                    f"ALT allele {self.alt!r} contains unsupported base {nucleotide!r}"
                )
        for nucleotide in self.ref:
            if nucleotide not in valid_nucleotides:
                raise InvalidAlleleError(
                    f"REF allele {self.ref!r} contains unsupported base {nucleotide!r}"
                )

    def __str__(self):
        return str({
            "chrom": self.chrom,
            "pos": self.pos,
            "ref": self.ref,
            "alt": self.alt,
            "b": self.b,
            "se": self.se,
            "nlog_pval": self.nlog_pval,
            "n": self.n,
            "alt_freq": self.alt_freq,
            "study_id": self.study_id,
            "ncase": self.ncase,
            "ncontrol": self.ncontrol,  # Added ncontrol
            "neff": self.neff,          # Added neff
            "imp_info": self.imp_info,
            "imp_z": self.imp_z,
            "vcf_filter": self.vcf_filter,
        })

    @staticmethod
    def read_from_file(
        input_file_path,
        fasta,
        chrom_col_num,
        pos_col_num,
        ea_col_num,
        nea_col_num,
        effect_col_num,
        se_col_num,
        pval_col_num,
        delimiter,
        header,
        ncase_col_num=None,
        study_id_col_num=None,
        ea_af_col_num=None,
        nea_af_col_num=None,
        imp_z_col_num=None,
        imp_info_col_num=None,
        ncontrol_col_num=None,
        alias=None,
    ):
        logging.info(f"Reading summary stats and mapping to FASTA: {input_file_path}")
        logging.debug(f"File path: {input_file_path}")
        logging.debug(f"CHR field: {chrom_col_num}")
        logging.debug(f"POS field: {pos_col_num}")
        logging.debug(f"EA field: {ea_col_num}")
        logging.debug(f"NEA field: {nea_col_num}")
        logging.debug(f"Effect field: {effect_col_num}")
        logging.debug(f"SE field: {se_col_num}")
        logging.debug(f"P fields: {pval_col_num}")
        logging.debug(f"Delimiter: {delimiter}")
        logging.debug(f"Header: {header}")
        logging.debug(f"ncase Field: {ncase_col_num}")
        logging.debug(f"Study variant identifier field: {study_id_col_num}")
        logging.debug(f"EA AF Field: {ea_af_col_num}")
        logging.debug(f"NEA AF Field: {nea_af_col_num}")
        logging.debug(f"IMP Z Score Field: {imp_z_col_num}")
        logging.debug(f"IMP INFO Field: {imp_info_col_num}")
        logging.debug(f"N Control Field: {ncontrol_col_num}")

        metadata = {
            "TotalVariants": 0,
            "VariantsNotRead": 0,
            "HarmonisedVariants": 0,
            "VariantsNotHarmonised": 0,
            "SwitchedAlleles": 0,
            "NormalisedVariants": 0,
        }
        file_idx = {}
        file_name, file_extension = os.path.splitext(input_file_path)

        if file_extension == ".gz":
            logging.info("Reading gzip file")
            f_handle = gzip.open(input_file_path, "rt")
        else:
            logging.info("Reading plain text file")
            f_handle = open(input_file_path)

        if header:
            logging.info(f"Skipping header: {f_handle.readline().strip()}")

        results = tempfile.TemporaryFile()
        p_value_handler = PvalueHandler()

        for line in f_handle:
            metadata["TotalVariants"] += 1
            columns = line.strip().split(delimiter)
            logging.debug(f"Input row: {columns}")

            try:
                if alias is not None:
                    if columns[chrom_col_num] in alias:
                        chrom = alias[columns[chrom_col_num]]
                    else:
                        chrom = columns[chrom_col_num]
                else:
                    chrom = columns[chrom_col_num]
            except Exception as exception_name:
                logging.debug(f"Skipping {columns}: {exception_name}")
                metadata["VariantsNotRead"] += 1
                continue

            try:
                pos = int(float(columns[pos_col_num]))
                if pos <= 0:
                    raise ValueError("Variant position must be a positive integer")
            except Exception as exception_name:
                logging.debug(f"Skipping {columns}: {exception_name}")
                metadata["VariantsNotRead"] += 1
                continue

            ref = str(columns[nea_col_num]).strip().upper()
            alt = str(columns[ea_col_num]).strip().upper()

            if ref == alt:
                logging.debug(f"Skipping: ref={ref} is the same as alt={alt}")
                metadata["VariantsNotRead"] += 1
                continue

            try:
                b = float(columns[effect_col_num])
            except Exception as exception_name:
                logging.debug(f"Skipping {columns}: {exception_name}")
                metadata["VariantsNotRead"] += 1
                continue

            try:
                se = float(columns[se_col_num])
            except Exception as exception_name:
                logging.debug(f"Skipping {columns}: {exception_name}")
                metadata["VariantsNotRead"] += 1
                continue

            try:
                pval = p_value_handler.parse_string(columns[pval_col_num])
                nlog_pval = p_value_handler.neg_log_of_decimal(pval)
            except Exception as exception_name:
                logging.debug(f"Skipping line {columns}, {exception_name}")
                metadata["VariantsNotRead"] += 1
                continue

            try:
                if ea_af_col_num is not None:
                    alt_freq = float(columns[ea_af_col_num])
                elif nea_af_col_num is not None:
                    alt_freq = 1 - float(columns[nea_af_col_num])
                else:
                    alt_freq = None
            except (IndexError, TypeError, ValueError) as exception_name:
                logging.debug(f"Could not parse allele frequency: {exception_name}")
                alt_freq = None

            try:
                study_id = columns[study_id_col_num].strip()
                if not study_id or study_id == ".":
                    study_id = None
            except (IndexError, TypeError, ValueError, AttributeError) as exception_name:
                logging.debug(f"Could not parse study variant identifier: {exception_name}")
                study_id = None

            try:
                ncase = float(columns[ncase_col_num]) if ncase_col_num is not None else None
            except (IndexError, TypeError, ValueError) as exception_name:
                logging.debug(f"Could not parse number of cases: {exception_name}")
                ncase = None

            try:
                ncontrol = float(columns[ncontrol_col_num]) if ncontrol_col_num is not None else None
            except (IndexError, TypeError, ValueError) as exception_name:
                logging.debug(f"Could not parse number of controls: {exception_name}")
                ncontrol = None

            try:
                n = (
                    ncase + ncontrol if ncase is not None and ncontrol is not None
                    else ncontrol if ncontrol is not None
                    else ncase
                )
            except (TypeError, ValueError) as exception_name:
                logging.debug(f"Could not determine total sample size: {exception_name}")
                n = None

            try:
                neff = (
                    4 / (1/ncase + 1/ncontrol) if ncase is not None and ncontrol is not None
                    else ncontrol if ncontrol is not None
                    else ncase
                )
            except (TypeError, ValueError, ZeroDivisionError) as exception_name:
                logging.debug(f"Could not calculate neff: {exception_name}")
                neff = None


            try:
                ncase = int(ncase) if ncase is not None else None
            except (TypeError, ValueError) as exception_name:
                logging.debug(f"Could not convert ncase to int: {exception_name}")
                ncase = None

            try:
                imp_info = float(columns[imp_info_col_num]) if imp_info_col_num is not None else None
            except (IndexError, TypeError, ValueError) as exception_name:
                logging.debug(f"Could not parse imputation INFO: {exception_name}")
                imp_info = None

            try:
                imp_z = float(columns[imp_z_col_num]) if imp_z_col_num is not None else None
            except (IndexError, TypeError, ValueError) as exception_name:
                logging.debug(f"Could not parse imputation Z score: {exception_name}")
                imp_z = None

            result = Gwas(
                chrom,
                pos,
                ref,
                alt,
                b,
                se,
                nlog_pval,
                n,
                alt_freq,
                study_id,
                ncase,
                ncontrol,
                neff,
                imp_info,
                imp_z,
            )

            logging.debug(f"Extracted row: {result}")

            try:
                result.check_alleles_are_valid()
            except InvalidAlleleError as exception_name:
                logging.debug(f"Skipping {columns}: {exception_name}")
                metadata["VariantsNotRead"] += 1
                continue

            try:
                result.check_reference_allele(fasta)
            except ReferenceAlleleMismatchError:
                try:
                    result.reverse_sign()
                    result.check_reference_allele(fasta)
                    metadata["SwitchedAlleles"] += 1
                except ReferenceAlleleMismatchError as exception_name:
                    logging.debug(f"Could not harmonise {columns}: {exception_name}")
                    metadata["VariantsNotHarmonised"] += 1
                    continue
            metadata["HarmonisedVariants"] += 1

            if len(ref) > 1 and len(alt) > 1:
                try:
                    result.normalise(fasta)
                except ReferenceSequenceFetchError:
                    raise
                except Exception as exception_name:
                    logging.debug(f"Could not normalise {columns}: {exception_name}")
                    metadata["VariantsNotHarmonised"] += 1
                    continue
                metadata["NormalisedVariants"] += 1

            if result.chrom not in file_idx:
                file_idx[result.chrom] = []
            heappush(file_idx[result.chrom], (result.pos, results.tell()))

            try:
                pickle.dump(result, results)
            except Exception as exception_name:
                logging.error(f"Could not write to {tempfile.gettempdir()}:", exception_name)
                raise exception_name

        f_handle.close()

        logging.info(f'Total variants: {metadata["TotalVariants"]}')
        logging.info(f'Variants could not be read: {metadata["VariantsNotRead"]}')
        logging.info(f'Variants harmonised: {metadata["HarmonisedVariants"]}')
        logging.info(
            f'Variants discarded during harmonisation: {metadata["VariantsNotHarmonised"]}'
        )
        logging.info(f'Alleles switched: {metadata["SwitchedAlleles"]}')
        logging.info(f'Normalised variants: {metadata["NormalisedVariants"]}')
        logging.info(
            f'Skipped {metadata["VariantsNotRead"] + metadata["VariantsNotHarmonised"]} of {metadata["TotalVariants"]}'
        )
        if (metadata["VariantsNotRead"] + metadata["VariantsNotHarmonised"]) / metadata[
            "TotalVariants"
        ] > 0.2:
            logging.warning(
                "More than 20% of variants not read or harmonised. Check your input"
            )

        return results, file_idx, metadata
