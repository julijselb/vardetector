from setuptools import setup, find_packages

VERSION = "0.0.1" 
DESCRIPTION = "Using bionumpy to detect DNA called somatic variants in RNA data."
LONG_DESCRIPTION = "Package that uses bionumpy to detect DNA called somatic variants (provided in a form of a .vcf file) in RNA data (.bam/.bai) and produces a report (1 variant per row) of the observed evidence."

# Setting up
setup(
       # the name must match the folder name 'verysimplemodule'
        name="vardetector", 
        version=VERSION,
        author="Julij Šelb",
        author_email="julij.selb@klinika-golnik.si",
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        packages=find_packages(),
        install_requires=["bionumpy"], # add any additional packages that 
        # needs to be installed along with your package. Eg: 'caer'

        keywords=["python", "vardetector", "somatic variants", "RNA data"],
        classifiers= [
            "Development Status :: 3 - Alpha",
            "Intended Audience :: Education",
            "Programming Language :: Python :: 2",
            "Programming Language :: Python :: 3",
            "Operating System :: MacOS :: MacOS X",
            "Operating System :: Microsoft :: Windows",
        ]
)