# Sparse Phenotyping Design
## ONE-CGIAR - BrIN - BSU – CIMMYT - Ángela Pacheco, Juan Burgueño
### 6th April 2022

#### Definition

In multi-environmental trials a sparse design is a phenotyping design in which not all the treatments are evaluated in all environments.

#### Description

The application generates a balanced sparse design. Genotypes are assigned to a given number of sets, there is a full set whose genotypes are evaluated in all environments while the remaining sets (sparse sets) are evaluated only in some environments. The sparse sets are of the same size and evaluated in the same number of environments. Checks are evaluated in all environments.
There are two options to allocate genotypes in sets given they are classified in separate groups, “Balanced sets” or “Balanced groups”. Using “Balanced sets” each set is compose of genotypes from all groups proportionally to group and set sizes, while “Balanced groups” allocate genotypes in the same group to the same set. Using the last options, comparison between genotypes in the same group will be larger than comparison between genotypes in separate groups. The full set always have genotypes from all groups following the size proportionality.
Table 1 shows a descriptive diagram of one sparse phenotyping design for 70 genotypes (including checks) tested in 5 environments. The 70 genotypes were separated in 6 sets, set 1 with 20 genotypes (including the checks) which are evaluated in all environments. Five sparse sets of the same size (10 genotypes) to be allocated in 4 environments each one. The 60 genotypes in each environment are evaluated in an alpha(0,1) experimental design. In future developments, other options of experimental design will be added.

#### Installation and start

After downloading the zip file, unzip it in a convenient place, a folder “\AppToolSparsePheno” will be created. To start the application, click on the run.vbs file. The application run into web browser program.

| Sets | Env1 | Env2 | Env3 | Env4 | Env5 |
| ---  |:----:|:----:|:----:|:----:|:----:|
| 1 (Full set) | `20 ` | `20 ` | `20 ` | `20 ` | `20 ` |
| 2  | `10 ` | `10 ` | `10 ` | `10 ` |  |
| 3  | `10 ` | `10 ` | `10 ` |  | `10 ` |
| 4  | `10 ` | `10 ` |  | `10 ` | `10 ` |
| 5  | `10 ` |  | `10 ` | `10 ` | `10 ` |
| 6  |  | `10 ` | `10 ` | `10 ` | `10 ` |

Table 1. Diagram of sparse phenotyping design for 70 genotypes and 5 environments. Numbers in cells are number of genotypes.


#### Data preparation

To use the application, it is necessary to have a data file in CSV format with at least three columns, GID, Role and Group. GID is a genotype identifier, role is a column with two values, Genotype and Check, and Groups is an identifier column of the group of genotypes. For check, the group label is “All”. Groups usually will be different families, different testers, or any other criteria for grouping genotypes. Columns names are mandatory for these three. More columns with additional information could be included and will be included in the final “Field Book Design”.

#### User guide

After starting the application, click on “Load file” box below “Input genotype names”. Automatically the program will produce sparse designs considering the number of genotypes, checks and groups in the file. The list of design is shown below “Summary options for available design”.
You can select a file name for the output design and the method to allocate genotypes to sets. The application is accompanied with some example files located in the folder “\AppToolSparsePheno\Examples”
The “optional parameters” section allows the user to filter the designs for the following parameters: Number of Sites, Number of Sites for un-replicated sets, Number of Replicates, Percent of saving, Plots by site, Block size, Number of columns in field, and Number of rows in field. The code 99999 generates designs for pre-specified values of the parameters. Regarding the total number of plots and the block size, the designs are arranged with different layouts in the field, i.e., number of rows and columns.
The application automatically adds checks to get a feasible number of genotypes in the experiment to build a balanced design with blocks of the same size. These checks could be substituted by genotypes, but information must be added manually into the design or generate a new list of genotypes with the new total number of genotypes to generate again the design.
Browsing the generated designs, you can select the one which fit better for your experimental conditions. Below “Options for get design” you can select the number of the desired design, click on “Get design” box and the design will be generated.
Two excel format files are saved in the same folder from which the list of genotypes was loaded. A file with the different options of design and the design selected and generated. You can generate more than one design, but this must be done one by one and taking care of changing the output file name.
The output design file has four sheets. Summary sheet has the values of the main parameters of the model. The efficiency compares the A-optimality of the alpha(0,1) design in each environment with the A-optimality of a complete block design. The “FieldBook” sheet has the experimental design for each environment in a tabular format and the “FieldMap” sheet has a diagram of the distribution of genotypes (numbers) and blocks (distinct colors) in the field.

To generate a sparse design, the following parameters must be considered, and some relationships between them must happen.

## Parameters
Sm: Number of environments

Km: Number of environments in which genotypes in sparse sets are evaluated.

H: Number of genotypes, excluding checks

C: Number of checks

N: total number of genotypes, including checks

R: Number of replicates of genotypes in each site

BbS: Number of blocks per site

B: Number of plots per block, block size

PbS: Total number of plots per environment

Q: Number of sets

noset: Number of sparse sets in each environment

Z1: Number of genotypes in the full set

Z2_q: Number of genotypes in the sparse sets 

## Relationships

Km < Sm

N = H + C

Q = Comb(Sm, Km) + 1, where Comb is for combinatorial

noset = Comb(Sm-1, Km-1)

(Q -1) * Z2_q + Z1 + C = H + C = N

PbS = (Z2_q * noset + Z1 + C)*R = B * BbS




