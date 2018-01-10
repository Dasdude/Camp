from parsers.dataset_make import heatmap,multiline_plot
import pandas as pd
import numpy as np
from parsers.df_data_extract import get_opposite_direction_dataset,get_same_direction_dataset,get_df_from_lane_diff
from parsers.plotmethodlibrary import *
import os.path  as path
import os
from matplotlib import pyplot as plt
if __name__ == '__main__':
	filename = 'alldone_density.csv'
	experiment_path = '../Dataset/DSRC/'
	plot_path  = path.join('..','Plots')
	if path.exists(plot_path)==-1: os.mkdir(plot_path)
	df_filepath = path.join(experiment_path, filename)
	df = pd.read_csv(df_filepath, index_col=False)
	df_same = get_same_direction_dataset(df)
	df_opposite = get_opposite_direction_dataset(df)
	densities_low = [0]
	# for idx,density_high in enumerate([70]):
	# 	tempdf_density = df_same[df_same.Av_Density >= densities_low[idx] ]
	# 	tempdf_density = tempdf_density[tempdf_density.Av_Density < density_high ]
	# 	tempdf_density.to_csv(path.join(experiment_path,'same_density_{}to{}.csv'.format(densities_low[idx],density_high)))
	for idx, density_high in enumerate([70]):
		tempdf_density = df_opposite[df_opposite.Av_Density >= densities_low[idx]]
		tempdf_density = tempdf_density[tempdf_density.Av_Density < density_high]
		tempdf_density.to_csv(path.join(experiment_path, 'opposite_density_{}to{}.csv'.format(densities_low[idx], density_high)))
		# plot_lane_difference_RSS_distance(df_same,path.join(plot_path,'lane diff density agg'),False,lane_set_cardinality=2)
	# plot_lane_difference_RSS_distance_density_sep(df=df_same, dir=path.join(plot_path, 'lane diff density_cardinality4'),show=False,
	#                                               lane_set_cardinality=4)
	# heatmap('', x_range=df['long'].values, y_theta=df['lat'].values, color_param=df['contact_lane'].values, show=True, save=False, dir=dir,
	#         xtitle='longitude', ytitle='latitude', colorbar_label='RSS(dbm)', cmap='tab20')
	cardinality=2;density_bins=3;
	mirror= True
	# for c in [2,4,8]:
	# 	for mirror in [True,False]:
	# 		plot_lane_difference_RSS_distance_density_sep_sender_ego(df=df_same, dir=path.join(plot_path,'RSS Dif-Final', 'Same Direction '
	# 		                                                                                                             'packets-lane '
	# 		                                                                                                            'signed diff '
	# 		                                                                                              'density_cardinality{}_Mirror:{'
	# 		                                                                                              '}_sender_ego Density_bins:{}'.format(c,
	# 		                                                                                                                              mirror,density_bins)),
	# 		                                                         show=False,
	# 		                                              lane_set_cardinality=c,mirror_lane_numbering=mirror,density_bins=density_bins)
	# density_step=25;window_ratio = 4
	# for lane_dif_card in [16]:
	# 	for c in [4]:
	# 		for mirror in [True,False]:
	# 			same_dir = path.join(plot_path, 'SlidingFinal', 'Same', 'Lane dif {}'.format(lane_dif_card), 'Mirror{}'.format(mirror),
	# 			                     'cardinality{}'.format(c), 'density sliding '
	# 			                                                'window_cardinality{' \
	# 			                                                '}_Mirror:{'
	# 			                                                '}_sender_ego_step:{},windowRatio{} '.format(c, mirror, density_step, window_ratio))
	# 			opp_dir = path.join(plot_path, 'SlidingFinal', 'Opposite', 'Lane dif {}'.format(lane_dif_card), 'Mirror{}'.format(mirror),
	# 			                     'cardinality{}'.format(c), 'density sliding '
	# 			                                                'window_cardinality{' \
	# 			                                                '}_Mirror:{'
	# 			                                                '}_sender_ego_step:{},windowRatio{} '.format(c, mirror, density_step, window_ratio))
	# 			# plot_lane_difference_RSS_distance_density_sep_sender_ego_window_density(df=df_same,
	# 			#                                                                         dir=same_dir, show=False,
	# 			#                                                                         lane_set_cardinality=c, mirror_lane_numbering=mirror,
	# 			#                                                                         density_bins=density_bins,density_step=density_step,
	# 			#                                                                         density_window_ratio=window_ratio,
	# 			#                                                                         lane_dif_group_card=lane_dif_card)
	# 			plot_lane_difference_RSS_distance_density_sep_sender_ego_window_density(df=df_opposite,
	# 			                                                                        dir=opp_dir, show=False,
	# 			                                                                        lane_set_cardinality=c, mirror_lane_numbering=mirror,
	# 			                                                                        density_bins=density_bins, density_step=density_step,
	# 			                                                                        density_window_ratio=window_ratio,
	# 			                                                                        lane_dif_group_card=lane_dif_card)
	# plt.pause(1000)