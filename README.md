# A Python package to detect DNA called variants (.vcf file) in RNA reads (.sam file)

# Installation
## install polars and scipy
python -m pip install polars
python -m pip install scipy

## clone the repository to the folder of your choice
/folder/of/your/choice git clone https://github.com/julijselb/vardetector.git

## export the path to vardetector for PYTHONPATH
export PYTHONPATH=/folder/of/your/choice/vardetector/vardetector

## you are ready to GO!

# Usage
put create_report_df import to the beginning of your py script

"
from detector import create_report_df

df: pl.DataFrame = create_report_df(/path/to/sam, /path/to/vcf)

"

# Acknowledgments and Funding

The development of the program was supported by the Slovenian Research and Innovation Agency (the vast majority (more than 50%) of the funding originated from grant J3-4516; other funding sources include grant J3-50109 and research programs P3-0360 and P1-0245).