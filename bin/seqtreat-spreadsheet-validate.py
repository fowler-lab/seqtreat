#! /usr/bin/python


import argparse
import pandas, re, numpy
from dateutil.parser import parse

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--spreadsheet",default=None,required=True,help="the path to the Seq&Treat spreadsheet to validate")
    parser.add_argument("--tables_path",default=None,required=True,help="the path to the Seq&Treat spreadsheet to validate")
    options = parser.parse_args()

    # load the spreadsheet
    df=pandas.read_excel(options.spreadsheet)

    print("%30s : %s" % ("spreadsheet",options.spreadsheet))
    print()

    SITES=pandas.read_csv(options.tables_path+"/SITES.csv")
    COUNTRIES=pandas.read_csv(options.tables_path+"/COUNTRIES_LOOKUP.csv")
    SEQUENCERS=pandas.read_csv(options.tables_path+"/SEQTREAT_SEQUENCERS.csv")
    AST_METHODS=pandas.read_csv(options.tables_path+'/SEQTREAT_AST_METHODS.csv')
    DRUG=pandas.read_csv(options.tables_path+"/DRUG_LOOKUP.csv.gz")

    # first work out which drugs are ok and which are not
    drug_list=[]
    failed_drugs=[]
    for i in df.columns:
        if 'method' in i:
            drug=i[:3]
            if drug.upper() not in DRUG.DRUG_ABBREVATION.unique():
                failed_drugs.append(drug)
            else:
                drug_list.append(drug)


    def validate_column(column_name,value):

        # the SITEID must exist in the table
        if column_name=='site_id':
            result=value in SITES.SITEID.unique()

        # as must the COUNTRY code
        elif column_name=='country_where_sample_taken':
            result=value in COUNTRIES.COUNTRY_CODE_3_LETTER.unique()

        elif column_name=='instrument_model':
            result=value in SEQUENCERS.instrument_model.unique()

        elif column_name=='drug_method':
            result=value in AST_METHODS.drug_method.unique()

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
            result=value in AST_METHODS.drug_method.unique()

        elif column_name=='phenotype':
            result=value in ['R','S','U']

        elif column_name=='cc':
            result=False
            if isinstance(value,str):
                result=value in ['WHO','UK']
            elif isinstance(value,float):
                result=value>0

        return result




    for column_name in ['dataset_name','site_id','subject_id','lab_id','isolate_number','sequence_replicate_number','country_where_sample_taken','collection_date','submission_date','ena_deposited','reads_file_1','reads_file_1_md5','reads_file_2','reads_file_2_md5','instrument_model','ena_run_accession','ena_sample_accession']:

        good_rows=0
        bad_rows=0
        bad_values=[]

        assert column_name in df.columns, column_name+" not found in spreadsheet!"

        for i in df[column_name]:

            if validate_column(column_name,i):
                good_rows+=1
            else:
                bad_rows+=1
                if i not in bad_values:
                    bad_values.append(i)

        if bad_rows==0:
            message="All rows pass validation"
        else:
            if len(bad_values)==1:
                message=str(good_rows)+ " rows which pass validation and "+str(bad_rows)+ " which fail validation. All failing values due to: "+str(bad_values[0])
            elif len(bad_values)==2:
                message=str(good_rows)+ " rows which pass validation and "+str(bad_rows)+ " which fail validation. All failing values due to: "+str(bad_values[0])+", "+str(bad_values[1])
            elif len(bad_values)==3:
                message=str(good_rows)+ " rows which pass validation and "+str(bad_rows)+ " which fail validation. All failing values due to: "+str(bad_values[0])+", "+str(bad_values[1])+", "+str(bad_values[2])
            else:
                message=str(good_rows)+ " rows which pass validation and "+str(bad_rows)+ " which fail validation. Number unique: "+str(len(bad_values))+ ". Example failing values: "+str(bad_values[0])+", "+str(bad_values[1])+", "+str(bad_values[2])

        print("%30s : %s" % (column_name,message))

    for drug_name in failed_drugs:
        print("%30s : %s" % (drug_name,"drug not recognised, please check"))


    for drug_name in drug_list:

        assert drug_name+"_method" in df.columns
        assert drug_name+"_cc" in df.columns
        assert drug_name+"_phenotype" in df.columns

        for field_name in ['method','phenotype','cc']:

            good_rows=0
            bad_rows=0
            bad_values=[]

            for i in df[drug_name+"_"+field_name]:

                if validate_column(field_name,i):
                    good_rows+=1
                else:
                    bad_rows+=1
                    if i not in bad_values:
                        bad_values.append(i)

            column_name=drug_name+"_"+field_name
            if bad_rows==0:
                message="All rows pass validation"
            else:
                if len(bad_values)==1:
                    message=str(good_rows)+ " rows which pass validation and "+str(bad_rows)+ " which fail validation. All failing values due to: "+str(bad_values[0])
                elif len(bad_values)==2:
                    message=str(good_rows)+ " rows which pass validation and "+str(bad_rows)+ " which fail validation. All failing values due to: "+str(bad_values[0])+", "+str(bad_values[1])
                elif len(bad_values)==3:
                    message=str(good_rows)+ " rows which pass validation and "+str(bad_rows)+ " which fail validation. All failing values due to: "+str(bad_values[0])+", "+str(bad_values[1])+", "+str(bad_values[2])
                else:
                    message=str(good_rows)+ " rows which pass validation and "+str(bad_rows)+ " which fail validation. Number unique: "+str(len(bad_values))+ ". Example failing values: "+str(bad_values[0])+", "+str(bad_values[1])+", "+str(bad_values[2])

            print("%30s : %s" % (column_name,message))
