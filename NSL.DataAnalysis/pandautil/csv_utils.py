import pandas as pd
import numpy as np

def get_pd_from_attr(data_frame,attr,value):
	# return dataframe where the attribute is equal to value
	return data_frame[getattr(data_frame,attr)==value]
