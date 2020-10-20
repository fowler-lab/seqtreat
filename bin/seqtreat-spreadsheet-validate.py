#! /usr/bin/python


import argparse, os

import pandas

import seqtreat

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--spreadsheet",default=None,required=True,help="the path to the Seq&Treat spreadsheet to validate")
    parser.add_argument("--tables_path",default=None,required=True,help="the path to the folder containing SITES.csv, COUNTRIES_LOOKUP.csv etc (cryptic-tables/ will work)")
    options = parser.parse_args()

    # load the spreadsheet
    donated_spreadsheet=pandas.read_excel(options.spreadsheet, dtype={'site_id':'object','collection_date':'object'})

    filestem=os.path.splitext(options.spreadsheet)[0]

    print("%30s : %s" % ("spreadsheet",options.spreadsheet))
    print()

    lookup_values={}

    SITES=pandas.read_csv(options.tables_path+"/SITES.csv")
    lookup_values['SITES']=SITES.SITEID.unique()

    COUNTRIES=pandas.read_csv(options.tables_path+"/COUNTRIES_LOOKUP.csv")
    lookup_values['COUNTRIES']=COUNTRIES.COUNTRY_CODE_3_LETTER.unique()

    SEQUENCERS=pandas.read_csv(options.tables_path+"/SEQTREAT_SEQUENCERS.csv")
    lookup_values['SEQUENCERS']=SEQUENCERS.instrument_model.unique()

    AST_METHODS=pandas.read_csv(options.tables_path+'/SEQTREAT_AST_METHODS.csv')
    lookup_values['AST_METHODS']=AST_METHODS.DRUG_METHOD.unique()

    DRUG=pandas.read_csv(options.tables_path+"/DRUG_CODES.csv")

    # first work out which drugs are ok and which are not
    drug_list=[]
    failed_drugs=[]
    for i in donated_spreadsheet.columns:
        if 'method' in i:
            drug=i[:3]
            if drug.upper() not in DRUG.DRUG_3_LETTER_CODE.unique():
                failed_drugs.append(drug)
            else:
                drug_list.append(drug)

    spreadsheet_pass=True

    for column_name in ['dataset_name','site_id','subject_id','lab_id','isolate_number','sequence_replicate_number','country_where_sample_taken','collection_date','submission_date','ena_deposited','reads_file_1','reads_file_1_md5','reads_file_2','reads_file_2_md5','instrument_model','ena_run_accession','ena_sample_accession']:

        good_rows=0
        bad_rows=0
        bad_values=[]

        assert column_name in donated_spreadsheet.columns, column_name+" not found in spreadsheet!"

        for i in donated_spreadsheet[column_name]:

            if seqtreat.validate_column(column_name,i,lookup_values):
                good_rows+=1
            else:
                bad_rows+=1
                if i not in bad_values:
                    bad_values.append(i)

        if bad_rows==0:
            message="All rows pass validation"
        else:
            spreadsheet_pass=False
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
        spreadsheet_pass=False

    for drug_name in drug_list:

        assert drug_name+"_method" in donated_spreadsheet.columns, drug_name+"_method not in columns"
        assert drug_name+"_cc" in donated_spreadsheet.columns, drug_name+"_cc not in columns"
        assert drug_name+"_phenotype" in donated_spreadsheet.columns, drug_name+"_phenotype not in columns"

        for field_name in ['method','phenotype','cc']:

            good_rows=0
            bad_rows=0
            bad_values=[]

            for i in donated_spreadsheet[drug_name+"_"+field_name]:

                if seqtreat.validate_column(field_name,i,lookup_values):
                    good_rows+=1
                else:
                    bad_rows+=1
                    if i not in bad_values:
                        bad_values.append(i)

            column_name=drug_name+"_"+field_name
            if bad_rows==0:
                message="All rows pass validation"
            else:
                spreadsheet_pass=False
                if len(bad_values)==1:
                    message=str(good_rows)+ " rows which pass validation and "+str(bad_rows)+ " which fail validation. All failing values due to: "+str(bad_values[0])
                elif len(bad_values)==2:
                    message=str(good_rows)+ " rows which pass validation and "+str(bad_rows)+ " which fail validation. All failing values due to: "+str(bad_values[0])+", "+str(bad_values[1])
                elif len(bad_values)==3:
                    message=str(good_rows)+ " rows which pass validation and "+str(bad_rows)+ " which fail validation. All failing values due to: "+str(bad_values[0])+", "+str(bad_values[1])+", "+str(bad_values[2])
                else:
                    message=str(good_rows)+ " rows which pass validation and "+str(bad_rows)+ " which fail validation. Number unique: "+str(len(bad_values))+ ". Example failing values: "+str(bad_values[0])+", "+str(bad_values[1])+", "+str(bad_values[2])

            print("%30s : %s" % (column_name,message))

    if not spreadsheet_pass:
        print("Spreadsheet DOES NOT PASS validation - please fix reported errors and repeat")
    else:
        print("Spreadsheet PASSES validation!")

        new_table_rows=[]

        for (idx,row) in donated_spreadsheet.iterrows():

            for drug in drug_list:

                if row[drug+"_phenotype"] in ['S','R','U']:

                    drug_name=drug.upper()

                    df=AST_METHODS[AST_METHODS.DRUG_METHOD==row[drug+"_method"]]

                    new_row={}

                    assert len(df)==1
                    new_row['UNIQUEID']='site.'+str(row['site_id'])+'.subj.'+str(row['subject_id'])+'.lab.'+str(row['lab_id'])+'.iso.'+str(row['isolate_number'])
                    new_row['DRUG'] = drug_name
                    new_row['METHOD_1']=df['METHOD_1'].values[0]
                    new_row['METHOD_2']=df['METHOD_2'].values[0]
                    new_row['METHOD_3']=df['METHOD_3'].values[0]
                    new_row['METHOD_CC']=row[drug+"_cc"]
                    new_row['METHOD_MIC']=''
                    new_row['PHENOTYPE']=row[drug+'_phenotype']
                    # df=df[['UNIQUEID','DRUG','METHOD_1','METHOD_2','METHOD_3','METHOD_CC','METHOD_MIC','PHENOTYPE']]
                    new_table_rows.append(new_row)

        df=pandas.DataFrame(new_table_rows)

        if not os.path.isfile(filestem+"_DST_MEASUREMENTS.csv"):
            df.to_csv(filestem+"_DST_MEASUREMENTS.csv")
        else:
            print("File already exists! Not saving...")
