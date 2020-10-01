#! /usr/bin/env python

import re, numpy
from dateutil.parser import parse

def validate_column(column_name,value,lookup_values):
    """Validates columns found in Seq&Treat tuberculosis AST donation spreadsheets.

    This function understands either the format of a passed column or uses values
    derived from lookup Pandas dataframes to check each value in a spreadsheet.

    Args:
        column_name (str): the name of the column. Not checked at present!
        value: the contents to check

    Returns:
        True/False
    """

    # the SITEID must exist in the table
    if column_name=='site_id':
        result=value in lookup_values['SITES']

    # as must the COUNTRY code
    elif column_name=='country_where_sample_taken':
        result=value in lookup_values['COUNTRIES']

    elif column_name=='instrument_model':
        result=value in lookup_values['SEQUENCERS']

    elif column_name=='isolate_number':
        try:
            result=value>0
        except:
            result=False

    elif column_name=='sequence_replicate_number':
        result=bool(re.match('^[_0-9]+$',str(value)))

    elif column_name in ['dataset_name','lab_id','subject_id']:
        result=bool(re.match('^[_\-A-Za-z0-9]+$',str(value)))

    elif column_name in ['collection_date','submission_date']:
        try:
            result=bool(parse(value))
        except:
            result=False

    elif column_name=='reads_file_1':
        result=bool(re.match('^[\-_A-Za-z0-9]+_R1.fastq.gz$',str(value)))

    elif column_name=='reads_file_2':
        result=bool(re.match('^[\-_A-Za-z0-9]+_R2.fastq.gz$',str(value)))

    elif column_name in ['reads_file_1_md5','reads_file_2_md5']:
        result=bool(re.match('^[a-z0-9]+$',str(value)))

    elif column_name in ['ena_deposited']:
        result=value in [True,False]

    elif column_name in ['ena_run_accession']:
        result=False
        if isinstance(value,float) and numpy.isnan(value):
            result=True
        elif isinstance(value,str):
            result=bool(re.match('^(E|D|S)RR[0-9]{6,}$',value))

    elif column_name in ['ena_sample_accession']:
        result=False
        if isinstance(value,float) and numpy.isnan(value):
            result=True
        elif isinstance(value,str):
            result=bool(re.match('^(E|D|S)RS[0-9]{6,}$',value))

    elif column_name=='method':
        result=value in lookup_values['AST_METHODS']

    elif column_name=='phenotype':
        result=value in ['R','S','U']

    elif column_name=='cc':
        result=False
        if isinstance(value,str):
            result=value in ['WHO','UK']
        elif isinstance(value,float):
            result=value>0

    return result
