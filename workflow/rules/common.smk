__author__ = "Jessika Nordin"
__copyright__ = "Copyright 2022, Martin Rippin"
__email__ = "jessika.nordin@scilifelab.uu.se"
__license__ = "GPL-3"

import pandas
import yaml

from hydra_genetics.utils.resources import load_resources
from hydra_genetics.utils.samples import *
from hydra_genetics.utils.units import *
from snakemake.utils import min_version
from snakemake.utils import validate

min_version("6.10.0")


### Set and validate config file
configfile: "config.yaml"
validate(config, schema="../schemas/config.schema.yaml")


### Load and validate resources file
config = load_resources(config, config["resources"])
validate(config, schema="../schemas/resources.schema.yaml")


### Read and validate samples file
samples = pandas.read_table(config["samples"], dtype=str).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")


### Read and validate units file
units = pandas.read_table(config["units"], dtype=str).set_index(["sample", "type", "flowcell", "lane"], drop=False).sort_index()
validate(units, schema="../schemas/units.schema.yaml")


### Set wildcard constraints
wildcard_constraints:
    sample="|".join(get_samples(samples)),
    unit="N|T|R",
    read="fastq[1|2]",
