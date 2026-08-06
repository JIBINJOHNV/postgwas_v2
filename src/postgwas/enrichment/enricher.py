# https://maayanlab.cloud/Enrichr/#libraries

import json
import requests
import time
import pandas as pd
from pathlib import Path

# --- CONFIGURATION ---
BASE_URL = 'https://maayanlab.cloud/Enrichr'

# Your exact list of libraries
ALL_LIBRARIES = [
    "Achilles_fitness_decrease", "Achilles_fitness_increase", "Aging_Perturbations_from_GEO_down",
    "Aging_Perturbations_from_GEO_up", "Allen_Brain_Atlas_10x_scRNA_2021", "Allen_Brain_Atlas_down",
    "Allen_Brain_Atlas_up", "ARCHS4_Cell-lines", "ARCHS4_IDG_Coexp", "ARCHS4_Kinases_Coexp",
    "ARCHS4_TFs_Coexp", "ARCHS4_Tissues", "Azimuth_2023", "Azimuth_Cell_Types_2021", "BioCarta_2013",
    "BioCarta_2015", "BioCarta_2016", "BioPlanet_2019", "BioPlex_2017", "Cancer_Cell_Line_Encyclopedia",
    "CCLE_Proteomics_2020", "CellMarker_2024", "CellMarker_Augmented_2021", "ChEA_2013", "ChEA_2015",
    "ChEA_2016", "ChEA_2022", "Chromosome_Location", "Chromosome_Location_hg19", "ClinVar_2019",
    "ClinVar_2025", "CM4AI_U2OS_Protein_Localization_Assemblies", "COMPARTMENTS_Curated_2025",
    "COMPARTMENTS_Experimental_2025", "CORUM", "COVID-19_Related_Gene_Sets", "COVID-19_Related_Gene_Sets_2021",
    "Data_Acquisition_Method_Most_Popular_Genes", "dbGaP", "DepMap_CRISPR_GeneDependency_CellLines_2023",
    "DepMap_WG_CRISPR_Screens_Broad_CellLines_2019", "DepMap_WG_CRISPR_Screens_Sanger_CellLines_2019",
    "Descartes_Cell_Types_and_Tissue_2021", "DGIdb_Drug_Targets_2024", "Diabetes_Perturbations_GEO_2022",
    "Disease_Perturbations_from_GEO_down", "Disease_Perturbations_from_GEO_up", "Disease_Signatures_from_GEO_down_2014",
    "Disease_Signatures_from_GEO_up_2014", "DisGeNET", "Drug_Perturbations_from_GEO_2014",
    "Drug_Perturbations_from_GEO_down", "Drug_Perturbations_from_GEO_up", "DrugMatrix", "DSigDB",
    "Elsevier_Pathway_Collection", "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X", "ENCODE_Histone_Modifications_2013",
    "ENCODE_Histone_Modifications_2015", "ENCODE_TF_ChIP-seq_2014", "ENCODE_TF_ChIP-seq_2015",
    "Enrichr_Libraries_Most_Popular_Genes", "Enrichr_Submissions_TF-Gene_Coocurrence", "Enrichr_Users_Contributed_Lists_2020",
    "Epigenomics_Roadmap_HM_ChIP-seq", "ESCAPE", "FANTOM6_lncRNA_KD_DEGs", "GeDiPNet_2023",
    "Gene_Perturbations_from_GEO_down", "Gene_Perturbations_from_GEO_up", "Genes_Associated_with_NIH_Grants",
    "GeneSigDB", "Genome_Browser_PWMs", "GlyGen_Glycosylated_Proteins_2022", "GO_Biological_Process_2021",
    "GO_Biological_Process_2023", "GO_Biological_Process_2025", "GO_Cellular_Component_2021",
    "GO_Cellular_Component_2023", "GO_Cellular_Component_2025", "GO_Molecular_Function_2021",
    "GO_Molecular_Function_2023", "GO_Molecular_Function_2025", "GTEx_Aging_Signatures_2021",
    "GTEx_Tissue_Expression_Down", "GTEx_Tissue_Expression_Up", "GTEx_Tissues_V8_2023", "GWAS_Catalog_2019",
    "GWAS_Catalog_2023", "GWAS_Catalog_2025", "HDSigDB_Human_2021", "HDSigDB_Mouse_2021", "HMDB_Metabolites",
    "HMS_LINCS_KinomeScan", "HomoloGene", "HuBMAP_ASCT_plus_B_augmented_w_RNAseq_Coexpression",
    "HuBMAP_ASCTplusB_augmented_2022", "Human_Gene_Atlas", "Human_Phenotype_Ontology", "HumanCyc_2015",
    "HumanCyc_2016", "huMAP", "IDG_Drug_Targets_2022", "InterPro_Domains_2019", "JASPAR_PWM_Human_2025",
    "JASPAR_PWM_Mouse_2025", "Jensen_COMPARTMENTS", "Jensen_DISEASES", "Jensen_DISEASES_Curated_2025",
    "Jensen_DISEASES_Experimental_2025", "Jensen_TISSUES", "KEA_2013", "KEA_2015", "KEGG_2013",
    "KEGG_2015", "KEGG_2016", "KEGG_2019_Human", "KEGG_2019_Mouse", "KEGG_2021_Human",
    "Kinase_Perturbations_from_GEO_down", "Kinase_Perturbations_from_GEO_up", "KOMP2_Mouse_Phenotypes_2022",
    "L1000_Kinase_and_GPCR_Perturbations_down", "L1000_Kinase_and_GPCR_Perturbations_up", "Ligand_Perturbations_from_GEO_down",
    "Ligand_Perturbations_from_GEO_up", "LINCS_L1000_Chem_Pert_Consensus_Sigs", "LINCS_L1000_Chem_Pert_down",
    "LINCS_L1000_Chem_Pert_up", "LINCS_L1000_CRISPR_KO_Consensus_Sigs", "LINCS_L1000_Ligand_Perturbations_down",
    "LINCS_L1000_Ligand_Perturbations_up", "lncHUB_lncRNA_Co-Expression", "MAGMA_Drugs_and_Diseases",
    "MAGNET_2023", "MCF7_Perturbations_from_GEO_down", "MCF7_Perturbations_from_GEO_up",
    "Metabolomics_Workbench_Metabolites_2022", "MGI_Mammalian_Phenotype_Level_4_2021",
    "MGI_Mammalian_Phenotype_Level_4_2024", "Microbe_Perturbations_from_GEO_down", "Microbe_Perturbations_from_GEO_up",
    "miRTarBase_2017", "MoTrPAC_2023", "Mouse_Gene_Atlas", "MSigDB_Computational", "MSigDB_Hallmark_2020",
    "MSigDB_Oncogenic_Signatures", "NCI-60_Cancer_Cell_Lines", "NCI-Nature_2016", "NIBR_DRUGseq_2025_down",
    "NIBR_DRUGseq_2025_up", "NURSA_Human_Endogenous_Complexome", "Old_CMAP_down", "Old_CMAP_up", "OMIM_Disease",
    "OMIM_Expanded", "Orphanet_Augmented_2021", "PanglaoDB_Augmented_2021", "Panther_2015", "Panther_2016",
    "PerturbAtlas", "PerturbAtlas_MouseGenePerturbationSigs", "PerturbSeq_ReplogleK562", "PerturbSeq_ReplogleRPE1",
    "Pfam_Domains_2019", "Pfam_InterPro_Domains", "PFOCR_Pathways_2023", "PhenGenI_Association_2021",
    "PheWeb_2019", "Phosphatase_Substrates_from_DEPOD", "PPI_Hub_Proteins", "Proteomics_Drug_Atlas_2023",
    "ProteomicsDB_2020", "Rare_Diseases_AutoRIF_ARCHS4_Predictions", "Rare_Diseases_AutoRIF_Gene_Lists",
    "Rare_Diseases_GeneRIF_ARCHS4_Predictions", "Rare_Diseases_GeneRIF_Gene_Lists", "Reactome_2022",
    "Reactome_Pathways_2024", "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO", "Rummagene_kinases",
    "Rummagene_signatures", "Rummagene_transcription_factors", "RummaGEO_DrugPerturbations_2025",
    "RummaGEO_GenePerturbations_2025", "Sciplex_Drug_Perturbation_Signatures_2025", "SILAC_Phosphoproteomics",
    "SubCell_BarCode", "SynGO_2022", "SynGO_2024", "SysMyo_Muscle_Gene_Sets", "Table_Mining_of_CRISPR_Studies",
    "Tabula_Muris", "Tabula_Sapiens", "TargetScan_microRNA", "TargetScan_microRNA_2017", "TF-LOF_Expression_from_GEO",
    "TF_Perturbations_Followed_by_Expression", "TG_GATES_2020", "The_Kinase_Library_2023", "The_Kinase_Library_2024",
    "Tissue_Protein_Expression_from_Human_Proteome_Map", "Tissue_Protein_Expression_from_ProteomicsDB",
    "TISSUES_Curated_2025", "TISSUES_Experimental_2025", "Transcription_Factor_PPIs", "TRANSFAC_and_JASPAR_PWMs",
    "TRRUST_Transcription_Factors_2019", "UK_Biobank_GWAS_v1", "Virus-Host_PPI_P-HIPSTer_2020",
    "Virus_Perturbations_from_GEO_down", "Virus_Perturbations_from_GEO_up", "VirusMINT", "WikiPathway_2021_Human",
    "WikiPathway_2023_Human", "WikiPathways_2013", "WikiPathways_2015", "WikiPathways_2016",
    "WikiPathways_2019_Human", "WikiPathways_2019_Mouse", "WikiPathways_2024_Human", "WikiPathways_2024_Mouse"
]


def enrichr_add_list(gene_list, description="Enrichr Analysis"):
    """
    Uploads a gene list to Enrichr and returns the userListId.
    """
    url = 'https://maayanlab.cloud/Enrichr/addList'
    genes_str = '\n'.join(gene_list)
    # Enrichr expects 'multipart/form-data' for the list
    payload = {
        'list': (None, genes_str),
        'description': (None, description)
    }
    try:
        response = requests.post(url, files=payload)
        response.raise_for_status()
        data = response.json()
        return data.get('userListId')
    except Exception as e:
        print(f"Error adding list: {e}")
        return None

def enrichr_get_results(user_list_id, library_name):
    """
    Fetches enrichment results for a specific library.
    """
    url = f'https://maayanlab.cloud/Enrichr/enrich?userListId={user_list_id}&backgroundType={library_name}'
    try:
        response = requests.get(url)
        response.raise_for_status()
        data = response.json()
        
        # The API returns a dictionary: { "LibraryName": [ [Rank, Term, Pval, ...], ... ] }
        if library_name in data:
            return data[library_name]
        else:
            print(f"Warning: No data returned for library {library_name}")
            return []
    except Exception as e:
        print(f"Error fetching results for {library_name}: {e}")
        return []



def run_enrichr(
    gene_list,
    libraries=ALL_LIBRARIES,
    output_dir=None,
    sample_id=None,
    save=True,
):
    """
    Run Enrichr analysis across multiple libraries.
    Args:
        gene_list (list): List of gene symbols.
        libraries (list): Enrichr libraries (e.g., ['GO_Biological_Process_2021'])
        output_dir (str | Path | None): Directory to save Enrichr results.
        sample_id (str | None): Sample/dataset identifier for filenames.
        save (bool): Whether to save output to disk.
    Returns:
        pd.DataFrame: Combined Enrichr results.
    """
    print(f"📤 Uploading {len(gene_list)} genes to Enrichr…")
    user_list_id = enrichr_add_list(gene_list)
    if not user_list_id:
        print("❌ Failed to upload gene list to Enrichr.")
        return pd.DataFrame()
    all_rows = []
    for lib in libraries:
        print(f"\t\t\t\t📊 Fetching results from: {lib}")
        results = enrichr_get_results(user_list_id, lib)
        # Enrichr Response Format:
        # 0: Rank
        # 1: Term Name
        # 2: P-value
        # 3: Z-score
        # 4: Combined Score
        # 5: Genes (list)
        # 6: Adjusted P-value (FDR)
        for item in results:
            all_rows.append({
                "Library": lib,
                "Rank": item[0],
                "Term": item[1],
                "PValue": item[2],
                "AdjPValue": item[6],
                "ZScore": item[3],
                "CombinedScore": item[4],
                "Genes": ", ".join(item[5]),
            })
        # Polite delay (important for Enrichr API)
        time.sleep(0.5)
    df = pd.DataFrame(all_rows)
    # ---------------------------------------------------------
    # Save output (optional)
    # ---------------------------------------------------------
    if save and output_dir is not None:
        df = df.sort_values(by=["Library", "AdjPValue"])
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        if sample_id is None:
            sample_id = "enrichr"
        out_file = output_dir / f"{sample_id}_enrichr_results.tsv"
        df.to_csv(out_file, sep="\t", index=False)
        print(f"💾 Enrichr results saved to: {out_file}")
    return df


# --- Main Execution ---
if __name__ == "__main__":
    # 1. Input Genes
    my_genes = ["TP53", "BRCA1", "EGFR", "TNF", "APOE", "IL6", "AKT1"]
    
    # 2. Select Libraries
    # Common libraries: 'GO_Biological_Process_2023', 'KEGG_2021_Human', 
    # 'WikiPathway_2023_Human', 'Human_Phenotype_Ontology', 'ClinVar_2019'
    # 3. Run Analysis
    df = run_enrichr(my_genes, ALL_LIBRARIES)

    if not df.empty:
        # 4. Sort and Display
        # Sort by P-Value (Raw) to see top hits
        df_sorted = df.sort_values(by=["Library", "P-Value"])
        
        print("\nTop 5 Results per Library (Sorted by Raw P-Value):")
        
        # Group by library and print head of each
        for lib in ALL_LIBRARIES:
            subset = df_sorted[df_sorted["Library"] == lib]
            if not subset.empty:
                print(f"\n--- {lib} ---")
                print(f"{'Term':<50} {'P-Value':<10} {'Adj P-Val':<10}")
                print("-" * 75)
                for _, row in subset.head(5).iterrows():
                    term = (row['Term'][:47] + '..') if len(row['Term']) > 47 else row['Term']
                    print(f"{term:<50} {row['P-Value']:.2e}   {row['Adj P-Value']:.2e}")

        # OPTIONAL: Save to CSV
        # df_sorted.to_csv("enrichr_results.csv", index=False)
    else:
        print("No results found.")
        
        