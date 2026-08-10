



conda create -n postgwas python=3.8
conda activate postgwas

pip install -r /Users/JJOHN41/Documents/developing_software/postgwas_underdevelopment/postgwas/src/gwas2vcf/requirements.txt
pip install wheel setuptools==59.8.0 pip==23.3.1
pip install git+https://github.com/bioinformed/vgraph@v1.4.0#egg=vgraph
python main.py -h



wget http://www.christianbenner.com/finemap_v1.4.2_MacOSX.tgz

mamba install pandas polars click  pyarrow scipy

git clone --recurse-submodules https://github.com/samtools/htslib.git
git clone https://github.com/samtools/bcftools.git


wget https://github.com/samtools/bcftools/releases/download/1.22/bcftools-1.22.tar.bz2
wget https://github.com/samtools/htslib/releases/download/1.22.1/htslib-1.22.1.tar.bz2

tar xjvf bcftools-1.22.tar.bz2
tar xjvf htslib-1.22.1.tar.bz2

rm bcftools-1.22.tar.bz2 htslib-1.22.1.tar.bz2

cd bcftools-1.22
wget -P plugins http://raw.githubusercontent.com/freeseek/score/master/{score.{c,h},{munge,liftover,metal,blup}.c,pgs.{c,mk}}

make
make plugins

cp bcftools ~/bin/
export BCFTOOLS_PLUGINS="/Users/JJOHN41/Documents/developing_software/postgwas_underdevelopment/postgwas/src/bcftools/bcftools-1.22/plugins"
export PATH="$HOME/bin:$PATH"

source ~/.bashrc


export PATH="/Users/JJOHN41/Documents/software_resources/softwares/magma_v1.10_mac:$PATH"


install.packages('optparse')


wget -P $HOME/GRCh37 http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz




install.packages('optparse')



mamba install bioconda::pyliftover
mamba install conda-forge::xgboost
mamba install conda-forge::fastparquet
pip install rich
pip install psutil





pip uninstall -y postgwas
rm -rf ~/.local/lib/python*/site-packages/postgwas*
rm -rf ~/miniconda3/envs/postgwas/lib/python*/site-packages/postgwas*
rm -rf postgwas.egg-info
pip install -e .



apt-get update && apt-get install -y cmake

git clone https://github.com/JianYang-Lab/GCTA.git
cd GCTA/
git submodule update --init
apt-get update && apt-get install -y build-essential ninja-build

docker build --platform=linux/amd64 -t jibinjv/postgwas:1.0 .
docker build --no-cache --platform=linux/amd64 -t jibinjv/postgwas:1.0 .

docker build --platform=linux/amd64 -t jibinjv/postgwas:1.4 .
docker build --no-cache --platform=linux/amd64 -t jibinjv/postgwas:1.4 .


docker build --platform=linux/amd64 -t jibinjv/postgwas:1.5 .
