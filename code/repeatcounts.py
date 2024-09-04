import sys,re
import pandas as pd
import numpy as np
from statsmodels.nonparametric.smoothers_lowess import lowess




def define_strata(stratum):#get stratum features
    feature={}
    if re.match("Gw.genome",stratum):
        feature["chr"] = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr13","chr14","chr15","chr16","chr17","chr18","chr20","chr21"]
        feature["group"] = "Gw.genome"
        feature["type"] = "Gw.genome"
        feature["length"] = 451874567

    elif re.match("Gw.S5.12",stratum):
        feature["chr"] = ["chr12"]
        feature["coordinates"] = ["5766000:11169000"]
        feature["group"] = "S5"
        feature["type"] = "chrX"
        feature["length"] = 5403000
    elif re.match("Gw.S5.Y",stratum):
        feature["chr"] = ["chrY"]
        feature["coordinates"] = ["2007000:7872000"]
        feature["group"] = "S5"
        feature["type"] = "chrY"
        feature["length"] = 5865000

    elif re.match("Gw.S4.2.12",stratum):
        feature["chr"] = ["chr12"]
        feature["coordinates"] = ["11169000:14597000"]
        feature["group"] = "S4.2"
        feature["type"] = "chrX"
        feature["length"] = 3428000
    elif re.match("Gw.S4.2.Y",stratum):
        feature["chr"] = ["chrY"]
        feature["coordinates"] = ["11728000:14720000"]
        feature["group"] = "S4.2"
        feature["type"] = "chrY"
        feature["length"] = 2992000

    elif re.match("Gw.Low.fst.12",stratum):
        feature["chr"] = ["chr12"]
        feature["coordinates"] = ["14597000:16615000"]
        feature["group"] = "Low.fst"
        feature["type"] = "chrX"
        feature["length"] = 2018000
    elif re.match("Gw.Low.fst.Y",stratum):
        feature["chr"] = ["chrY"]
        feature["coordinates"] = ["10240000:11728000"]
        feature["group"] = "Low.fst"
        feature["type"] = "chrY"
        feature["length"] = 1488000

    elif re.match("Gw.S4.1.12",stratum):
        feature["chr"] = ["chr12"]
        feature["coordinates"] = ["16615000:20650000"]
        feature["group"] = "S4.1"
        feature["type"] = "chrX"
        feature["length"] = 4035000
    elif re.match("Gw.S4.1.Y",stratum):
        feature["chr"] = ["chrY"]
        feature["coordinates"] = ["7872000:10240000","200000:2007000","15371000:15613000"]
        feature["group"] = "S4.1"
        feature["type"] = "chrY"
        feature["length"] = 4417000

    elif re.match("Gw.S3.12",stratum):
        feature["chr"] = ["chr12"]
        feature["coordinates"] = ["20650000:23704000"]
        feature["group"] = "S3"
        feature["type"] = "chrX"
        feature["length"] = 3054000
    elif re.match("Gw.S3.Y",stratum):
        feature["chr"] = ["chrY"]
        feature["coordinates"] = ["18806000:22713000"]
        feature["group"] = "S3"
        feature["type"] = "chrY"
        feature["length"] = 3907000




    elif re.match("Ga.genome",stratum):
        feature["chr"] = ["NC_053212.1","NC_053213.1","NC_053214.1","NC_053215.1","NC_053216.1","NC_053217.1","NC_053218.1","NC_053219.1","NC_053220.1","NC_053221.1","NC_053222.1","NC_053224.1","NC_053225.1","NC_053226.1","NC_053227.1","NC_053228.1","NC_053229.1","NC_053231.1","NC_053232.1"]
        feature["group"] = "Ga.genome"
        feature["type"] = "Ga.genome"
        feature["length"] = 395156579

    elif re.match("Ga.S5",stratum):
        feature["chr"] = ["NC_053223.1"]
        feature["coordinates"] = ["4400000:9650000"]
        feature["group"] = "S5"
        feature["type"] = "Ga.chr12"
        feature["length"] = 5250000

    elif re.match("Ga.S4.2",stratum):
        feature["chr"] = ["NC_053223.1"]
        feature["coordinates"] = ["9650000:13120000"]
        feature["group"] = "S4.2"
        feature["type"] = "Ga.chr12"
        feature["length"] = 3470000

    elif re.match("Ga.Low.fst",stratum):
        feature["chr"] = ["NC_053223.1"]
        feature["coordinates"] = ["13120000:14414000","14919000:15180000"]
        feature["group"] = "Low.fst"
        feature["type"] = "Ga.chr12"
        feature["length"] = 1555000

    elif re.match("Ga.S4.1",stratum):
        feature["chr"] = ["NC_053223.1"]
        feature["coordinates"] = ["14414000:14919000","15180000:18643000"]
        feature["group"] = "S4.1"
        feature["type"] = "Ga.chr12"
        feature["length"] = 3968000

    elif re.match("Ga.S3",stratum):
        feature["chr"] = ["NC_053223.1"]
        feature["coordinates"] = ["18643000:20700000"]
        feature["group"] = "S3"
        feature["type"] = "Ga.chr12"
        feature["length"] = 2007000

    return feature


def extract_strata(df,stratum_feature): # extract data in a given region
    

    if re.search("genome",stratum_feature["type"]):
        chr_list = stratum_feature["chr"]
        group = stratum_feature["group"]
        type = stratum_feature["type"]
        length = stratum_feature["length"]

        df_chr_filtered = df[df['Chromosome'].isin(chr_list)]
        data = pd.DataFrame(columns=['ID', 'Chromosome', 'Start', 'End', 'RepeatType', 'K2P_Distance', 'DIV_Distance']) #Create a empty dataframe
        data = pd.concat([data,df_chr_filtered],ignore_index=True)
        data['group']=group
        data['type']=type
        data['length']=length

    else:
        chr_list = stratum_feature["chr"]
        coordinates = stratum_feature["coordinates"]
        group = stratum_feature["group"]
        type = stratum_feature["type"]
        length = stratum_feature["length"]

        df_chr_filtered = df[df['Chromosome'].isin(chr_list)]
        data = pd.DataFrame(columns=['ID', 'Chromosome', 'Start', 'End', 'RepeatType', 'K2P_Distance', 'DIV_Distance']) #Create a empty dataframe
        
        for pos in coordinates:
            start = pos.split(":")[0]
            end = pos.split(":")[1]
            extracted=df_chr_filtered.loc[(df_chr_filtered['Start'] > int(start)) & (df_chr_filtered['End'] < int(end))]
            data=pd.concat([data,extracted],ignore_index=True)

        
        data['group']=group
        data['type']=type
        data['length']=length

    return data     


def normalize_and_smooth(bin_counts_data,stratum_length):
    bin_counts_data.loc[:,'Normalized_Count']= bin_counts_data['Count'] / stratum_length * 1000000
    smoothed_values = lowess(bin_counts_data['Normalized_Count'], bin_counts_data['Midpoint'], frac=0.2)
    bin_counts_data['Smoothed_Normalized_Count'] = smoothed_values[:, 1]

    return bin_counts_data




def K2P_bin_summarize(data):

    # Define bins from 0 to 60 with a step of 2
    bins = np.arange(-1, 69, 2)
    stratum = data.loc[0,['group']]['group']
    sex = data.loc[0,['type']]['type']


    # Use pd.cut to categorize K2P_Distance into bins and count occurrences in each bin
    data['K2P_Binned'] = pd.cut(data['K2P_Distance'], bins=bins, right=False)

    # Count the number of occurrences in each bin
    bin_counts = data['K2P_Binned'].value_counts().sort_index()

    # Convert results to a DataFrame
    bin_counts_data = bin_counts.reset_index()
    bin_counts_data.columns = ['K2P_Range', 'Count']

    # Calculate the midpoint for each bin range
    bin_counts_data['Midpoint'] = bin_counts_data['K2P_Range'].apply(lambda x: (x.left + x.right) / 2)

    # Reorder columns to place 'Midpoint' between 'K2P_Range' and 'Count'
    bin_counts_data = bin_counts_data[['K2P_Range', 'Midpoint', 'Count']]

    stratum_length = data.loc[0,['length']]['length']

    result=normalize_and_smooth(bin_counts_data,stratum_length)
    result['Type']=sex
    result['Grouop']=stratum

    return result


def main():
    if len(sys.argv) != 3:
        print("Usage: python k2p_bin_count.py <input_file_path> <output_file_path>")
        sys.exit(1)
    input=sys.argv[1]
    output=sys.argv[2]

    try:
        # Read the file into a DataFrame
        df = pd.read_csv(input, sep='\t', header=None)
    except Exception as e:
        print(f"Error reading file: {e}")
        sys.exit(1)

    df.columns = ['ID', 'Chromosome', 'Start', 'End', 'RepeatType', 'K2P_Distance', 'DIV_Distance'] # Add header   

    strata_Gw = ["Gw.S5.12",
            "Gw.S5.Y",
            "Gw.S4.2.12",
            "Gw.S4.2.Y",
            "Gw.Low.fst.12",
            "Gw.Low.fst.Y",
            "Gw.S4.1.12",
            "Gw.S4.1.Y",
            "Gw.S3.12",
            "Gw.S3.Y",
            "Gw.genome"
    ]

    strata_Ga = ["Ga.S5",
            "Ga.S4.2",
            "Ga.Low.fst",
            "Ga.S4.1",
            "Ga.S3",
            "Ga.genome"
    ]   


    strata=strata_Ga

    final=pd.DataFrame(columns=["K2P_Range","Midpoint","Count","Normalized_Count","Smoothed_Normalized_Count","type"])
    
    for stratum in strata:
        print(f'{stratum} start:')
        stratum_feature = define_strata(stratum)
        data = extract_strata(df,stratum_feature)
        result = K2P_bin_summarize(data)
        final=pd.concat([final,result],ignore_index=True)
        print(f'{stratum} done!')

    try:
        final.to_csv(output, sep='\t', index=False)
        print(f"Results have been written to {output}")
    except Exception as e:
        print(f"Error writing to file: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()



