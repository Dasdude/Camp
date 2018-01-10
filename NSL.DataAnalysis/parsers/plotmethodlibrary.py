import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import os.path as path
import matplotlib.colors as colors
import matplotlib.cm as cmx


def multiline_plot_color_bar(title, x_range_list, y_theta_list, line_id_list, xtitle='', ytitle='', dir='', show=True, save=True, fun_handle=np.mean,
                      bins=80,
                   new_fig=True, line_indicator='--', y_range_bound=[-100, 30], err_bar=True, legend_label_prefix='Lane:'):
	def running_stat(X, Y, total_bins, fun_handle):
		total_bins
		bins = np.linspace(X.min(), X.max(), total_bins)
		bins = np.arange(0, 800, 10)
		delta = bins[1] - bins[0]
		centers = []
		idx = np.digitize(X, bins)
		running_median = []
		variance = []
		for k in np.arange(total_bins):
			val_points = Y[idx == k]
			if len(val_points) == 0:
				continue
			else:
				running_median += [fun_handle(Y[idx == k])]
				variance += [np.std(Y[idx == k])]
				centers += [bins[k]]
		return centers - delta / 2, running_median, variance

	if new_fig:
		fig = plt.figure(figsize=(15, 10))
		ax2 = fig.add_subplot(111)
		plt.title(title)
		jet = cm = plt.get_cmap('jet')
		cNorm = colors.Normalize(vmin=0, vmax=50)
		scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
		scalarMap.set_array([])
	count = 0
	for x, y, id in zip(x_range_list, y_theta_list, line_id_list):
		if len(x) == 0 or len(y) == 0:
			continue
		x_bin, y_med, y_var = running_stat(x, y, bins, fun_handle)
		# plt.plot(x_bin, y_med, line_indicator, label='Lane:{}'.format(id))
		if title.find('Samples') == -1:
			plt.ylim(-100, -30)
			plt.xlim(0, 800)
		else:
			y_var = None
		if err_bar:
			plt.errorbar(x_bin, y_med, yerr=y_var, label='{}{}'.format(legend_label_prefix, id))
		else:
			colorVal = scalarMap.to_rgba(id)
			plt.plot(x_bin, y_med, line_indicator, label='{}:{}'.format(legend_label_prefix, id), color=colorVal)
		# plt.plot(x_bin, y_med, line_indicator, label='{}:{}'.format(legend_label_prefix, id))
	plt.colorbar(scalarMap)
	plt.hold(True)
	# try:
	# 	plt.legend()
	# except:
	# 	print('error')
	plt.xlabel(xtitle)
	plt.ylabel(ytitle)
	if (not dir == '') and save:
		plt.savefig(os.path.join(dir, title + '.png'))
	plt.show(block=False) if show else plt.close()


def multiline_plot(title, x_range_list, y_theta_list, line_id_list, xtitle='', ytitle='', dir='', show=True, save=True, fun_handle=np.mean,bins=80,
                   new_fig=True,line_indicator='--',y_range_bound=[-100,30],err_bar=True,legend_label_prefix='Lane:'):
	def running_stat(X, Y, total_bins, fun_handle):
		total_bins
		bins = np.linspace(X.min(), X.max(), total_bins)
		bins = np.arange(0,800,10)
		delta = bins[1] - bins[0]
		centers = []
		idx = np.digitize(X, bins)
		running_median = []
		variance = []
		for k in np.arange(total_bins):
			val_points = Y[idx == k]
			if len(val_points) == 0:
				continue
			else:
				running_median += [fun_handle(Y[idx == k])]
				variance += [np.std(Y[idx == k])]
				centers += [bins[k]]
		return centers - delta / 2, running_median,variance
	if new_fig:
		fig = plt.figure(figsize=(15, 10))
		ax2 = fig.add_subplot(111)
		plt.title(title)
		jet = cm = plt.get_cmap('jet')
		cNorm = colors.Normalize(vmin=0, vmax=40)
		scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
		scalarMap.set_array([])
	count = 0
	for x, y, id in zip(x_range_list, y_theta_list, line_id_list):
		if len(x) == 0 or len(y) == 0:
			continue
		x_bin, y_med,y_var = running_stat(x, y, bins, fun_handle)
		# plt.plot(x_bin, y_med, line_indicator, label='Lane:{}'.format(id))
		if title.find('Samples')==-1:
			plt.ylim(-100,-30)
			plt.xlim(0, 800)
		else:
			y_var=None
		if err_bar:
			plt.errorbar(x_bin, y_med,yerr=y_var, label='{}{}'.format(legend_label_prefix,id))
		else:
			# colorVal = scalarMap.to_rgba(id)
			# plt.plot(x_bin, y_med, line_indicator, label='{}:{}'.format(legend_label_prefix,id),color=colorVal)
			plt.plot(x_bin, y_med, line_indicator, label='{}:{}'.format(legend_label_prefix, id))
	# plt.colorbar(scalarMap)
	# plt.hold(True)
	try:
		plt.legend()
	except:
		print('error')
	plt.xlabel(xtitle)
	plt.ylabel(ytitle)
	if (not dir == '') and save:
		plt.savefig(os.path.join(dir, title + '.png'))
	plt.show(block=False) if show else plt.close()


def plot__RSS_distance_lane_group_multiplot_density_seperate(df, dir, show, lane_set_cardinality=1, mirror_lane_numbering=False, density_bins=5):
	def agg_plot(df_list, plot_group_name, dir):
		if not os.path.exists(dir):
			os.mkdir(dir)
		for i in ['Min', 'Max', 'Mean', 'Median', 'Samples']:
			if not path.exists(path.join(dir, i)):
				os.mkdir(path.join(dir, i))
		theta = []
		range = []
		RSS = []
		density = []
		ego = []
		for df in df_list:
			theta += [df['angle'].values]
			range += [df['distance'].values]
			RSS += [df['RSS'].values]
			density += [df['Total_density'].values]
			ego += [set(df['ego_lane'].values)]
		# Polar all_range
		multiline_plot('{}RSS-Distance Mean'.format(plot_group_name), x_range_list=range, y_theta_list=RSS, show=show, dir=path.join(dir, 'Mean'),
		               xtitle='Range(meter)', ytitle='RSS', fun_handle=np.mean, line_id_list=ego)
		multiline_plot('{}RSS-Distance Median'.format(plot_group_name), x_range_list=range, y_theta_list=RSS, show=show, dir=path.join(dir,
		                                                                                                                               'Median'),
		               xtitle='Range(meter)', ytitle='RSS', fun_handle=np.median, line_id_list=ego)
		multiline_plot('{}RSS-Distance Max'.format(plot_group_name), x_range_list=range, y_theta_list=RSS, show=show, dir=path.join(dir, 'Max'),
		               xtitle='Range(meter)', ytitle='RSS', fun_handle=np.max, line_id_list=ego)
		multiline_plot('{}RSS-Distance Min'.format(plot_group_name), x_range_list=range, y_theta_list=RSS, show=show, dir=path.join(dir, 'Min'),
		               xtitle='Range(meter)', ytitle='RSS', fun_handle=np.min, line_id_list=ego)
		multiline_plot('{}RSS-Distance Samples'.format(plot_group_name), x_range_list=range, y_theta_list=RSS, show=show,
		               dir=path.join(dir, 'Samples'), xtitle='Range(meter)', ytitle='Samples', fun_handle=len, line_id_list=ego)

	# multiline_plot('{} Density AGG RSS-Distance std'.format(plot_group_name), x_range_list=range, y_theta_list=RSS, show=show, dir=path.join(dir,
	#
	# 'Mean'),
	#                xtitle='Range(meter)', ytitle='RSS', fun_handle=np.std, line_id_list=ego)
	def plotter(df, plot_group_name, dir):
		if not os.path.exists(dir):
			os.mkdir(dir)
		theta = df['angle'].values
		range = df['distance'].values
		RSS = df['RSS'].values
		density = df['Total_density'].values
		# Polar all_range
		heatmap('{} Density Heat Map RSS-Distance'.format(plot_group_name), x_range=range, y_theta=RSS, color_param=density, proj='', show=show,
		        dir=dir, xtitle='Range(meter)', ytitle='RSS', colorbar_label='density(dbm)')

	if not os.path.exists(dir):
		os.mkdir(dir)

	max_density = df.Total_density.max()
	min_density = df.Total_density.min()
	step_density = (max_density - min_density) / density_bins
	density_list = min_density + (step_density * np.arange(density_bins + 1))
	for i in np.arange(density_bins):
		diff_df = []
		for e_lane in np.arange(16 / lane_set_cardinality):
			lane_dir = os.path.join(dir, str(e_lane))
			if mirror_lane_numbering:
				df_east = df[df.ego_lane < 8]
				df_west = df[df.ego_lane >= 8]

				tempdf_ego = pd.concat([df_east[(df_east.ego_lane // lane_set_cardinality) == e_lane],
				                        df_west[((15 - df_west.ego_lane) // lane_set_cardinality) == e_lane]])
			else:
				tempdf_ego = df[(df.ego_lane // lane_set_cardinality) == e_lane]
			tempdf_ego = tempdf_ego[tempdf_ego.Total_density < density_list[i + 1]]
			tempdf_ego = tempdf_ego[tempdf_ego.Total_density > density_list[i]]
			diff_df += [tempdf_ego]
		# plotter(tempdf_ego, 'Lane#{} Lane Difference:{}'.format(e_lane, dif), os.path.join(dir, 'Lane Diff:' + str(dif)))
		agg_plot(diff_df, 'DensityRange {}:{} '.format(int(density_list[i] / 16), int(density_list[i + 1] / 16)), os.path.join(dir))

def heatmap(title, x_range, y_theta, color_param, cmap='plasma', proj=None, xtitle='', ytitle='', dir='', show=True, save=True, s=.1,
            colorbar_label=''):
	fig = plt.figure(figsize=(15, 10))
	if proj == 'polar':
		ax2 = fig.add_subplot(111, projection='polar')
		plt.title(title)
		plt.scatter(y_theta, x_range, c=color_param, cmap=cmap, s=s)
	else:
		ax2 = fig.add_subplot(111)
		plt.title(title)
		plt.scatter(x_range, y_theta, c=color_param, cmap=cmap, s=s)

	cb = plt.colorbar()
	cb.set_label(colorbar_label)
	plt.xlabel(xtitle)
	plt.ylabel(ytitle)
	if (not dir == '') and save:
		plt.savefig(os.path.join(dir, title + '.png'))
	plt.show(block=False) if show else plt.close()

	return


def plot_lane_difference_RSS_distance(df, dir, show, lane_set_cardinality=1, mirror_lane_numbering=False,density_bins=5,plot_name=''):
	def agg_plot(df_list, plot_group_name, dir):
		if not os.path.exists(dir):
			os.mkdir(dir)
		for i in ['Min', 'Max', 'Mean', 'Median', 'Samples']:
			if not path.exists(path.join(dir, i)):
				os.mkdir(path.join(dir, i))
		theta = []
		range = []
		RSS = []
		density = []
		ego = []
		for df in df_list:
			theta += [df['angle'].values]
			range += [df['distance'].values]
			RSS += [df['RSS'].values]
			density += [df['Total_density'].values]
			ego += [set(df['ego_lane'].values)]
		# Polar all_range
		multiline_plot('{} RSS-Distance Mean'.format(plot_group_name), x_range_list=range, y_theta_list=RSS, show=show, dir=path.join(dir, 'Mean'),
		               xtitle='Range(meter)', ytitle='RSS', fun_handle=np.mean, line_id_list=ego)
		multiline_plot('{} RSS-Distance Median'.format(plot_group_name), x_range_list=range, y_theta_list=RSS, show=show,
		               dir=path.join(dir, 'Median'), xtitle='Range(meter)', ytitle='RSS', fun_handle=np.median, line_id_list=ego)
		multiline_plot('{} RSS-Distance Max'.format(plot_group_name), x_range_list=range, y_theta_list=RSS, show=show, dir=path.join(dir, 'Max'),
		               xtitle='Range(meter)', ytitle='RSS', fun_handle=np.max, line_id_list=ego)
		multiline_plot('{} RSS-Distance Min'.format(plot_group_name), x_range_list=range, y_theta_list=RSS, show=show, dir=path.join(dir, 'Min'),
		               xtitle='Range(meter)', ytitle='RSS', fun_handle=np.min, line_id_list=ego)
		multiline_plot('{} RSS-Distance Samples'.format(plot_group_name), x_range_list=range, y_theta_list=RSS, show=show,
		               dir=path.join(dir, 'Samples'), xtitle='Range(meter)', ytitle='Samples', fun_handle=len, line_id_list=ego)

	# multiline_plot('{} Density AGG RSS-Distance std'.format(plot_group_name), x_range_list=range, y_theta_list=RSS, show=show, dir=path.join(dir,
	#
	# 'Mean'),
	#                xtitle='Range(meter)', ytitle='RSS', fun_handle=np.std, line_id_list=ego)
	def plotter(df, plot_group_name, dir):
		if not os.path.exists(dir):
			os.mkdir(dir)
		theta = df['angle'].values
		range = df['distance'].values
		RSS = df['RSS'].values
		density = df['Total_density'].values
		# Polar all_range
		heatmap('{} Density Heat Map RSS-Distance'.format(plot_group_name), x_range=range, y_theta=RSS, color_param=density, proj='', show=show,
		        dir=dir, xtitle='Range(meter)', ytitle='RSS', colorbar_label='density(dbm)')

	if not os.path.exists(dir):
		os.mkdir(dir)
	max_density = df.Total_density.max()
	min_density = df.Total_density.min()
	step_density = (max_density - min_density) / density_bins
	density_list = min_density + (step_density * np.arange(density_bins + 1))
	for dif in np.arange(16):
		diff_df = []
		for e_lane in np.arange(16 / lane_set_cardinality):
			# lane_dir = os.path.join(dir, str(e_lane))
			# if not os.path.exists(lane_dir):
			# 	os.mkdir(lane_dir)
			if mirror_lane_numbering:
				df_east = df[df.ego_lane < 8]
				df_west = df[df.ego_lane >= 8]

				tempdf_ego = pd.concat([df_east[(df_east.ego_lane // lane_set_cardinality) == e_lane],
				                        df_west[((15 - df_west.ego_lane) // lane_set_cardinality) == e_lane]])
			else:
				tempdf_ego = df[(df.ego_lane // lane_set_cardinality) == e_lane]
			tempdf_ego = tempdf_ego[np.abs(tempdf_ego.ego_lane - tempdf_ego.contact_lane) == dif]
			diff_df += [tempdf_ego]
			# plotter(tempdf_ego, 'Lane#{} Lane Difference:{}'.format(e_lane, dif), os.path.join(dir, 'Lane Diff:' + str(dif)))
		agg_plot(diff_df, plot_name+'Lane Difference: {}'.format(dif), os.path.join(dir))

def get_seprate_based_on_lane_dif(df,mirror_lane_numbering,lane_set_cardinality):
	all_dif = 16*[[]]
	for dif in np.arange(16):
		for e_lane in np.arange(16 / lane_set_cardinality):
			# lane_dir = os.path.join(dir, str(e_lane))
			# if not os.path.exists(lane_dir):
			# 	os.mkdir(lane_dir)
			if mirror_lane_numbering:
				df_east = df[df.ego_lane < 8]
				df_west = df[df.ego_lane >= 8]

				tempdf_ego = pd.concat([df_east[(df_east.contact_lane // lane_set_cardinality) == e_lane],
				                        df_west[((15 - df_west.contact_lane) // lane_set_cardinality) == e_lane]])
			else:
				tempdf_ego = df[(df.contact_lane // lane_set_cardinality) == e_lane]
			tempdf_ego = tempdf_ego[np.abs(tempdf_ego.ego_lane - tempdf_ego.contact_lane) == dif]
			all_dif[dif] =all_dif[dif]+ [tempdf_ego]
	return all_dif


def plot_sliding_window_density(df, dir, show, density_steps=100,density_window_ratio=.1, plot_name=''):
	def agg_plot(df_list, plot_group_name, dir):
		if not os.path.exists(dir):
			os.makedirs(dir)
		for i in ['Min', 'Max', 'Mean', 'Median', 'Samples']:
			if not path.exists(path.join(dir, i)):
				os.makedirs(path.join(dir, i))
		theta = []
		range = []
		RSS = []
		density = []
		ego = []
		for df in df_list:
			theta += [df['angle'].values]
			range += [df['distance'].values]
			RSS += [df['RSS'].values]
			density += [df['Av_Density'].min()]
			ego += [set(df['contact_lane'].values)]
		# Polar all_range
		multiline_plot_color_bar('{} RSS-Distance Mean'.format(plot_group_name), x_range_list=range, y_theta_list=RSS, show=show, dir=path.join(dir, 'Mean'),
		               xtitle='Range(meter)', ytitle='RSS', fun_handle=np.mean, line_id_list=density,err_bar=False,legend_label_prefix='Density:')
		multiline_plot_color_bar('{} RSS-Distance Median'.format(plot_group_name), x_range_list=range, y_theta_list=RSS, show=show,
		               dir=path.join(dir, 'Median'), xtitle='Range(meter)', ytitle='RSS', fun_handle=np.median, line_id_list=density, err_bar=False,
		               legend_label_prefix='Density:')
		multiline_plot_color_bar('{} RSS-Distance Max'.format(plot_group_name), x_range_list=range, y_theta_list=RSS, show=show, dir=path.join(dir, 'Max'),
		               xtitle='Range(meter)', ytitle='RSS', fun_handle=np.max, line_id_list=density, err_bar=False, legend_label_prefix='Density:')
		multiline_plot_color_bar('{} RSS-Distance Min'.format(plot_group_name), x_range_list=range, y_theta_list=RSS, show=show, dir=path.join(dir, 'Min'),
		               xtitle='Range(meter)', ytitle='RSS', fun_handle=np.min, line_id_list=density, err_bar=False, legend_label_prefix='Density:')
		multiline_plot_color_bar('{} RSS-Distance Samples'.format(plot_group_name), x_range_list=range, y_theta_list=RSS, show=show,
		               dir=path.join(dir, 'Samples'), xtitle='Range(meter)', ytitle='Samples', fun_handle=len, line_id_list=density, err_bar=False,
		               legend_label_prefix='Density:')

	window_size = density_window_ratio
	max_density = df.Av_Density.max()-window_size
	min_density = df.Av_Density.min()+window_size
	# window_size = (max_density-min_density)*density_window_ratio
	step_density = (max_density - min_density) / density_steps
	density_list = min_density + (step_density * np.arange(density_steps + 1))
	df_list = []
	for density in density_list:
		tempdf_density = df[df.Av_Density >= density-window_size]
		tempdf_density = tempdf_density[tempdf_density.Av_Density < density+window_size]
		df_list+=[tempdf_density]
		# plotter(tempdf_ego, 'Lane#{} Lane Difference:{}'.format(e_lane, dif), os.path.join(dir, 'Lane Diff:' + str(dif)))
	agg_plot(df_list, plot_name + 'density', os.path.join(dir))
	print('hist in dir',dir,'called')
	plt.figure()
	plt.hist(df['Av_Density'],bins=50,range=(0,50))
	plt.legend()
	plt.savefig(os.path.join(dir, plot_name+'Density_Hist.png'))


def plot_lane_difference_RSS_distance_sender_ego(df, dir, show, lane_set_cardinality=1, mirror_lane_numbering=False, density_bins=5, plot_name=''):
	def agg_plot(df_list, plot_group_name, dir):
		if not os.path.exists(dir):
			os.mkdir(dir)
		for i in ['Min', 'Max', 'Mean', 'Median', 'Samples']:
			if not path.exists(path.join(dir, i)):
				os.mkdir(path.join(dir, i))
		theta = []
		range = []
		RSS = []
		density = []
		ego = []
		for df in df_list:
			theta += [df['angle'].values]
			range += [df['distance'].values]
			RSS += [df['RSS'].values]
			density += [df['Av_Density'].values]
			ego += [set(df['contact_lane'].values)]
		# Polar all_range
		multiline_plot('{} RSS-Distance Mean'.format(plot_group_name), x_range_list=range, y_theta_list=RSS, show=show, dir=path.join(dir, 'Mean'),
		               xtitle='Range(meter)', ytitle='RSS', fun_handle=np.mean, line_id_list=ego)
		multiline_plot('{} RSS-Distance Median'.format(plot_group_name), x_range_list=range, y_theta_list=RSS, show=show,
		               dir=path.join(dir, 'Median'), xtitle='Range(meter)', ytitle='RSS', fun_handle=np.median, line_id_list=ego)
		multiline_plot('{} RSS-Distance Max'.format(plot_group_name), x_range_list=range, y_theta_list=RSS, show=show, dir=path.join(dir, 'Max'),
		               xtitle='Range(meter)', ytitle='RSS', fun_handle=np.max, line_id_list=ego)
		multiline_plot('{} RSS-Distance Min'.format(plot_group_name), x_range_list=range, y_theta_list=RSS, show=show, dir=path.join(dir, 'Min'),
		               xtitle='Range(meter)', ytitle='RSS', fun_handle=np.min, line_id_list=ego)
		multiline_plot('{} RSS-Distance Samples'.format(plot_group_name), x_range_list=range, y_theta_list=RSS, show=show,
		               dir=path.join(dir, 'Samples'), xtitle='Range(meter)', ytitle='Samples', fun_handle=len, line_id_list=ego)
	if not os.path.exists(dir):
		os.mkdir(dir)
	max_density = df.Av_Density.max()
	min_density = df.Av_Density.min()
	step_density = (max_density - min_density) / density_bins
	density_list = min_density + (step_density * np.arange(density_bins + 1))

	for dif in np.arange(16):
		diff_df = []
		for e_lane in np.arange(16 / lane_set_cardinality):
			# lane_dir = os.path.join(dir, str(e_lane))
			# if not os.path.exists(lane_dir):
			# 	os.mkdir(lane_dir)
			if mirror_lane_numbering:
				df_east = df[df.ego_lane < 8]
				df_west = df[df.ego_lane >= 8]

				tempdf_ego = pd.concat([df_east[(df_east.contact_lane // lane_set_cardinality) == e_lane],
				                        df_west[((15 - df_west.contact_lane) // lane_set_cardinality) == e_lane]])
			else:
				tempdf_ego = df[(df.contact_lane // lane_set_cardinality) == e_lane]
			tempdf_ego = tempdf_ego[np.abs(tempdf_ego.ego_lane - tempdf_ego.contact_lane) == dif]
			diff_df += [tempdf_ego]
		# plotter(tempdf_ego, 'Lane#{} Lane Difference:{}'.format(e_lane, dif), os.path.join(dir, 'Lane Diff:' + str(dif)))
		print('LaneDiff',dif)
		agg_plot(diff_df, plot_name + 'Lane Difference: {}'.format(dif), os.path.join(dir))


def plot_RSS_distance_moving_density(df, dir, show, density_bins=100,lane_set_cardinality=16, plot_name='',bin_window=10):
	def agg_plot(df_list, plot_group_name, dir):
		if not os.path.exists(dir):
			os.mkdir(dir)
		for i in ['Min', 'Max', 'Mean', 'Median', 'Samples']:
			if not path.exists(path.join(dir, i)):
				os.mkdir(path.join(dir, i))
		theta = []
		range = []
		RSS = []
		density = []
		ego = []
		for df in df_list:
			theta += [df['angle'].values]
			range += [df['distance'].values]
			RSS += [df['RSS'].values]
			density += [df['Av_Density'].values]
			ego += [set(df['contact_lane'].values)]
		# Polar all_range
		multiline_plot('{} RSS-Distance Mean'.format(plot_group_name), x_range_list=range, y_theta_list=RSS, show=show, dir=path.join(dir, 'Mean'),
		               xtitle='Range(meter)', ytitle='RSS', fun_handle=np.mean, line_id_list=ego)
		multiline_plot('{} RSS-Distance Median'.format(plot_group_name), x_range_list=range, y_theta_list=RSS, show=show,
		               dir=path.join(dir, 'Median'), xtitle='Range(meter)', ytitle='RSS', fun_handle=np.median, line_id_list=ego)
		multiline_plot('{} RSS-Distance Max'.format(plot_group_name), x_range_list=range, y_theta_list=RSS, show=show, dir=path.join(dir, 'Max'),
		               xtitle='Range(meter)', ytitle='RSS', fun_handle=np.max, line_id_list=ego)
		multiline_plot('{} RSS-Distance Min'.format(plot_group_name), x_range_list=range, y_theta_list=RSS, show=show, dir=path.join(dir, 'Min'),
		               xtitle='Range(meter)', ytitle='RSS', fun_handle=np.min, line_id_list=ego)
		multiline_plot('{} RSS-Distance Samples'.format(plot_group_name), x_range_list=range, y_theta_list=RSS, show=show,
		               dir=path.join(dir, 'Samples'), xtitle='Range(meter)', ytitle='Samples', fun_handle=len, line_id_list=ego)

	# multiline_plot('{} Density AGG RSS-Distance std'.format(plot_group_name), x_range_list=range, y_theta_list=RSS, show=show, dir=path.join(dir,
	#
	# 'Mean'),
	#                xtitle='Range(meter)', ytitle='RSS', fun_handle=np.std, line_id_list=ego)
	def plotter(df, plot_group_name, dir):
		if not os.path.exists(dir):
			os.mkdir(dir)
		theta = df['angle'].values
		range = df['distance'].values
		RSS = df['RSS'].values
		density = df['Total_density'].values
		# Polar all_range
		heatmap('{} Density Heat Map RSS-Distance'.format(plot_group_name), x_range=range, y_theta=RSS, color_param=density, proj='', show=show,
		        dir=dir, xtitle='Range(meter)', ytitle='RSS', colorbar_label='density(dbm)')

	if not os.path.exists(dir):
		os.mkdir(dir)
	max_density = df.Av_Density.max()
	min_density = df.Av_Density.min()
	step_density = (max_density - min_density) / density_bins
	density_list = min_density + (step_density * np.arange(density_bins + 1))
	for idx,density in enumerate(density_list):
		density_df=[]
		density_df+=df[(df.Av_Density>density)and(df.Av_Density < density_list[np.minimum(idx+bin_window,len(density_list)-1)])]
		# plotter(tempdf_ego, 'Lane#{} Lane Difference:{}'.format(e_lane, dif), os.path.join(dir, 'Lane Diff:' + str(dif)))
		agg_plot(density_df, plot_name + 'Moving Density Window'.format(dif), os.path.join(dir))
def plot_signed_lane_difference_RSS_distance_sender_ego(df, dir, show, lane_set_cardinality=1, mirror_lane_numbering=False,
                                                                   density_bins=5,
                                                     plot_name=''):
	def agg_plot(df_list_pos,df_list_neg, plot_group_name, dir):
		if not os.path.exists(dir):
			os.mkdir(dir)
		for i in ['Min', 'Max', 'Mean', 'Median', 'Samples']:
			if not path.exists(path.join(dir, i)):
				os.mkdir(path.join(dir, i))
		theta = []
		range = []
		RSS = []
		density = []
		ego = []
		theta_neg = []
		range_neg = []
		RSS_neg = []
		density_neg = []
		ego_neg = []
		for df in df_list_pos:
			theta += [df['angle'].values]
			range += [df['distance'].values]
			RSS += [df['RSS'].values]
			density += [df['Av_Density'].values]
			ego += [set(df['contact_lane'].values)]
		for df in df_list_neg:
			theta_neg += [df['angle'].values]
			range_neg += [df['distance'].values]
			RSS_neg += [df['RSS'].values]
			density_neg += [df['Av_Density'].values]
			ego_neg += [set(df['contact_lane'].values)]
		# Polar all_range
		multiline_plot('{} RSS-Distance Mean'.format(plot_group_name), x_range_list=range, y_theta_list=RSS, show=show, dir=path.join(dir, 'Mean'),
		               xtitle='Range(meter)', ytitle='RSS', fun_handle=np.mean, line_id_list=ego,bins=75,new_fig=False,line_indicator='*')
		multiline_plot('{} RSS-Distance Mean'.format(plot_group_name), x_range_list=range_neg, y_theta_list=RSS_neg, show=show, dir=path.join(dir,
		                                                                                                                                  'Mean'),
		               xtitle='Range(meter)', ytitle='RSS', fun_handle=np.mean, line_id_list=ego_neg, bins=75, new_fig=True, line_indicator='--')
		# multiline_plot('{} RSS-Distance Median'.format(plot_group_name), x_range_list=range, y_theta_list=RSS, show=show,
		#                dir=path.join(dir, 'Median'), xtitle='Range(meter)', ytitle='RSS', fun_handle=np.median, line_id_list=ego)
		# multiline_plot('{} RSS-Distance Max'.format(plot_group_name), x_range_list=range, y_theta_list=RSS, show=show, dir=path.join(dir, 'Max'),
		#                xtitle='Range(meter)', ytitle='RSS', fun_handle=np.max, line_id_list=ego)
		# multiline_plot('{} RSS-Distance Min'.format(plot_group_name), x_range_list=range, y_theta_list=RSS, show=show, dir=path.join(dir, 'Min'),
		#                xtitle='Range(meter)', ytitle='RSS', fun_handle=np.min, line_id_list=ego)
		# multiline_plot('{} RSS-Distance Samples'.format(plot_group_name), x_range_list=range, y_theta_list=RSS, show=show,
		#                dir=path.join(dir, 'Samples'), xtitle='Range(meter)', ytitle='Samples', fun_handle=len, line_id_list=ego)

	# multiline_plot('{} Density AGG RSS-Distance std'.format(plot_group_name), x_range_list=range, y_theta_list=RSS, show=show, dir=path.join(dir,
	#
	# 'Mean'),
	#                xtitle='Range(meter)', ytitle='RSS', fun_handle=np.std, line_id_list=ego)
	def plotter(df, plot_group_name, dir):
		if not os.path.exists(dir):
			os.mkdir(dir)
		theta = df['angle'].values
		range = df['distance'].values
		RSS = df['RSS'].values
		density = df['Av_Density'].values
		# Polar all_range
		heatmap('{} Density Heat Map RSS-Distance'.format(plot_group_name), x_range=range, y_theta=RSS, color_param=density, proj='', show=show,
		        dir=dir, xtitle='Range(meter)', ytitle='RSS', colorbar_label='density(dbm)')

	if not os.path.exists(dir):
		os.mkdir(dir)
	max_density = df.Total_Density.max()
	min_density = df.Total_Density.min()
	step_density = (max_density - min_density) / density_bins
	density_list = min_density + (step_density * np.arange(density_bins + 1))
	for dif in np.arange(7):
		diff_df = []
		diff_df_neg = []
		for e_lane in np.arange(16 / lane_set_cardinality):
			# lane_dir = os.path.join(dir, str(e_lane))
			# if not os.path.exists(lane_dir):
			# 	os.mkdir(lane_dir)
			if mirror_lane_numbering:
				df_east = df[df.ego_lane < 8]
				df_west = df[df.ego_lane >= 8]

				tempdf_ego = pd.concat([df_east[(df_east.contact_lane // lane_set_cardinality) == e_lane],
				                        df_west[((15 - df_west.contact_lane) // lane_set_cardinality) == e_lane]])
			else:
				tempdf_ego = df[(df.contact_lane // lane_set_cardinality) == e_lane]
			tempdf_ego_neg = tempdf_ego[tempdf_ego.ego_lane - tempdf_ego.contact_lane == -dif]
			tempdf_ego = tempdf_ego[tempdf_ego.ego_lane - tempdf_ego.contact_lane == dif]
			diff_df += [tempdf_ego]
			diff_df_neg += [tempdf_ego]
		# plotter(tempdf_ego, 'Lane#{} Lane Difference:{}'.format(e_lane, dif), os.path.join(dir, 'Lane Diff:' + str(dif)))
		agg_plot(diff_df,diff_df_neg, plot_name + 'Lane Difference(reciever - transmitter): {}'.format(dif), os.path.join(dir))
def	plot_lane_difference_RSS_distance_density_sep(df, dir, show, lane_set_cardinality=1, mirror_lane_numbering=False, density_bins=5, plot_name=''):
	max_density = df.Av_Density.max()
	min_density = df.Av_Density.min()
	step_density = (max_density - min_density) / density_bins
	density_list = min_density + (step_density * np.arange(density_bins + 1))
	for binidx in np.arange(density_bins):
		tempdf_density = df[df.Av_Density < density_list[binidx + 1]]
		tempdf_density = tempdf_density[df.Av_Density > density_list[binidx]]
		os.makedirs(path.join(dir, 'density_' + str(int(density_list[binidx]))),exist_ok=True)
		plot_lane_difference_RSS_distance(tempdf_density,path.join(dir,'density_'+str(int(density_list[binidx]))),show=False,
		                                  lane_set_cardinality=lane_set_cardinality)

def plot_lane_difference_RSS_distance_density_sep_sender_ego(df, dir, show, lane_set_cardinality=1, mirror_lane_numbering=False, density_bins=5,
                                                  plot_name=''):
	max_density = df.Av_Density.max()
	min_density = df.Av_Density.min()
	step_density = (max_density - min_density) / density_bins
	density_list = min_density + (step_density * np.arange(density_bins + 1))
	for binidx in np.arange(density_bins):
		tempdf_density = df[df.Av_Density < density_list[binidx + 1]]
		tempdf_density = tempdf_density[tempdf_density.Av_Density > density_list[binidx]]
		os.makedirs(path.join(dir, 'density_' + str(int(density_list[binidx]))), exist_ok=True)
		plot_lane_difference_RSS_distance_sender_ego(tempdf_density, path.join(dir, 'density_' + str(int(density_list[binidx]))), show=False,
		                                  lane_set_cardinality=lane_set_cardinality,mirror_lane_numbering=mirror_lane_numbering)


def plot_lane_difference_RSS_distance_density_sep_sender_ego_window_density(df, dir, show, lane_set_cardinality=1, mirror_lane_numbering=False,
                                                                    density_bins=5,
                                                             plot_name='',density_step=100,density_window_ratio=.1,lane_dif_group_card=1):
	# df_list = get_seprate_based_on_lane_dif(df,mirror_lane_numbering,lane_set_cardinality)
	print('ello')
	for dif in np.arange(16//lane_dif_group_card):
		print(dif)
		for e_lane in np.arange(16 / lane_set_cardinality):
			print(e_lane)
			# lane_dir = os.path.join(dir, str(e_lane))
			# if not os.path.exists(lane_dir):
			# 	os.mkdir(lane_dir)
			if mirror_lane_numbering:
				df_east = df[df.ego_lane < 8]
				df_west = df[df.ego_lane >= 8]

				tempdf_ego = pd.concat([df_east[(df_east.contact_lane // lane_set_cardinality) == e_lane],
				                        df_west[((15 - df_west.contact_lane) // lane_set_cardinality) == e_lane]])
			else:
				tempdf_ego = df[(df.contact_lane // lane_set_cardinality) == e_lane]
			tempdf_ego = tempdf_ego[np.abs(tempdf_ego.ego_lane - tempdf_ego.contact_lane) >= dif*lane_dif_group_card]
			tempdf_ego = tempdf_ego[np.abs(tempdf_ego.ego_lane - tempdf_ego.contact_lane) < int((dif+1)*lane_dif_group_card)]
			# os.makedirs(path.join(dir, 'lane_dif' + str(int(dif)),str(int(e_lane))), exist_ok=True)
			plot_sliding_window_density(tempdf_ego, path.join(dir, 'lane_dif' + str(int(dif*lane_dif_group_card))+'to'+ str(
				int((dif+1) * lane_dif_group_card))),show=False,plot_name=str(set(
				tempdf_ego.contact_lane.values)),density_steps=density_step,density_window_ratio=density_window_ratio)
			print(dif)


def plot_lane_difference_RSS_distance_density_sep_sender_ego_sliding(df, dir, show, lane_set_cardinality=1, mirror_lane_numbering=False,
                                                                    density_bins=5,
                                                             plot_name=''):
	max_density = df.Av_Density.max()
	min_density = df.Av_Density.min()
	step_density = (max_density - min_density) / density_bins
	density_list = min_density + (step_density * np.arange(density_bins + 1))
	for binidx in np.arange(density_bins):
		tempdf_density = df[df.Av_Density < density_list[binidx + 1]]
		tempdf_density = tempdf_density[tempdf_density.Av_Density > density_list[binidx]]
		os.makedirs(path.join(dir, 'density_' + str(int(density_list[binidx]))), exist_ok=True)
		plot_lane_difference_RSS_distance_sender_ego(tempdf_density, path.join(dir, 'density_' + str(int(density_list[binidx]))), show=False,
		                                             lane_set_cardinality=lane_set_cardinality)
def plot_RSS_distance_moving_density_window(df):
	max_density = df.Av_Density.max()
	min_density = df.Av_Density.min()
	step_density = (max_density - min_density) / density_bins
	density_list = min_density + (step_density * np.arange(density_bins + 1))
	for binidx in np.arange(density_bins):
		tempdf_density = df[df.Av_Density < density_list[binidx + 1]]
		tempdf_density = tempdf_density[tempdf_density.Av_Density > density_list[binidx]]
		os.makedirs(path.join(dir, 'density_' + str(int(density_list[binidx]))), exist_ok=True)
		plot_lane_difference_RSS_distance_sender_ego(tempdf_density, path.join(dir, 'density_' + str(int(density_list[binidx]))), show=False,
		                                             lane_set_cardinality=lane_set_cardinality)