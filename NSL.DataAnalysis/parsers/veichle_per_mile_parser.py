import pandas as pd
import numpy as np
import pandautil.csv_utils as csv_u
FILE_NAME_TEMPLATE = '{}_Vehicles_Per_Mile_{}Bound.csv'
if __name__ == '__main__':
    path = '../Dataset/sf/' + FILE_NAME_TEMPLATE.format('GG', 'East')
    df = pd.read_csv(path)
    df.plot(y=['1', '2', '3', '4', '5', '6', '7', '8'], x='HM')

    print(df)
def split_col(data_frame:pd.DataFrame,attr:str,attr1,attr2,splitter):
    df[attr1] = df[attr].str.rpartition(splitter)[0]
    df[attr2] = df[attr].str.rpartition(splitter)[2]
    return df