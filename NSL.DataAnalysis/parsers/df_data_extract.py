import pandas as pd
import numpy as np
def get_same_direction_dataset(samples):
	samples_1 = samples[samples.ego_lane < 8]
	samples_2 = samples[samples.ego_lane > 7]
	samples_1_1 = samples_1[samples_1.contact_lane < 8]
	# samples_1_1['Av_Density']=samples_1_1['Av_East_Density']
	samples_1_1.loc[:, ('Av_Density')] = samples_1_1.loc[:, ('Av_East_Density')]
	# res = samples_1_1.loc[:, ('2')]
	# for i in np.arange(3,8):
	# 	res = res+samples_1_1.loc[:, (str(i))]
	# samples_1_1.loc[:, ('Av_Density')]=res/6
	samples_2_2 = samples_2[samples_2.contact_lane > 7]
	# samples_2_2['Av_Density'] = samples_2_2['Av_West_Density']
	samples_2_2.loc[:, ('Av_Density')] = samples_2_2.loc[:, ('Av_West_Density')]
	res = pd.concat([samples_1_1, samples_2_2], axis=0)
	print(res['Av_Density'].max())
	return res
def get_opposite_direction_dataset(samples):
	samples_1 = samples[samples.ego_lane < 8]
	samples_2 = samples[samples.ego_lane > 7]
	samples_1_2 = samples_1[samples_1.contact_lane > 7]
	# samples_1_2['Av_Density']=(samples_1_2['Av_East_Density']+ samples_1_2['Av_West_Density'])/2
	samples_1_2.loc[:,('Av_Density')]= (samples_1_2.loc[:,('Av_East_Density')] + samples_1_2.loc[:,('Av_West_Density')])/2
	samples_2_1 = samples_2[samples_2.contact_lane < 8]
	# samples_2_1['Av_Density'] = (samples_2_1['Av_East_Density'] + samples_2_1['Av_West_Density']) / 2
	samples_2_1.loc[:, ('Av_Density')] = (samples_2_1.loc[:, ('Av_East_Density')] + samples_2_1.loc[:, ('Av_West_Density')])/2
	res = pd.concat([samples_1_2, samples_2_1], axis=0)
	return res

def get_df_from_lane_diff(df,lane_dif):
	return df[np.abs(df.contact_lane-df.ego_lane)==lane_dif]

