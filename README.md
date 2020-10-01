# Validation of Seq&Treat data

Python code to validate the spreadsheets containing *M. tuberculosis* AST data donated to the project. Likely to expand in future.

## Installation

To install as a module so `seqtreat-spreadsheet-validate.py` is in yoru ``$PATH` (and will dynamically update if you do a $git pull$ to pick up future code changes), do

```
$ git clone https://github.com/philipwfowler/seqtreat.git
$ cd seqtreat
$ python setup.py develop --user
```

Then this should work

```
$ seqtreat-spreadsheet-validate.py --help
usage: seqtreat-spreadsheet-validate.py [-h] --spreadsheet SPREADSHEET
                                        --tables_path TABLES_PATH

optional arguments:
  -h, --help            show this help message and exit
  --spreadsheet SPREADSHEET
                        the path to the Seq&Treat spreadsheet to validate
  --tables_path TABLES_PATH
                        the path to the folder containing SITES.csv,
                        COUNTRIES_LOOKUP.csv etc (cryptic-tables/ will work)
```
