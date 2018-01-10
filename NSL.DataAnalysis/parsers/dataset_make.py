from parsers import DSCR
from multiprocessing import Pool
import os
import pandas as pd
import re
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
import sklearn.cluster as cl
import os.path as path
# cpid is the Id of owner of packet. so for RXE it is the senders ID same as uniqueobe id but 0 based
def sliding(df):
	fig, ax = plt.subplots()
	plt.subplots_adjust(left=0.25, bottom=0.25)
	density_step = 5
	angle_step=np.pi/10
	a0 = .4
	f0 = 30
	tdf = df[df.Total_density<f0]
	tdf = tdf[f0-30<tdf.Total_density ]
	tdf = tdf[tdf.angle < a0]
	tdf = tdf[a0-.3<tdf.angle]
	l = plt.scatter(tdf.distance.values, tdf.RSS.values, color='red',s=4)
	# plt.axis([0, 1, -10, 10])

	axcolor = 'lightgoldenrodyellow'
	axfreq = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
	axamp = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)

	sfreq = Slider(axfreq, 'density', 0, 500.0 ,valinit=f0)
	samp = Slider(axamp, 'angle', -np.pi, np.pi,valinit=a0)

	def update(val):
		amp = samp.val
		freq = sfreq.val
		tdf = df[df.Total_density < freq]
		tdf = tdf[freq - density_step < tdf.Total_density]
		tdf = tdf[tdf.angle < amp]
		tdf = tdf[amp - angle_step < tdf.angle]

		l.set_offsets(np.transpose(np.array([tdf.distance.values,tdf.RSS.values])))
		print(np.array([tdf.distance.values, tdf.RSS.values]).shape)
		fig.canvas.draw_idle()


	sfreq.on_changed(update)
	samp.on_changed(update)

	resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
	button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')

	def reset(event):
		sfreq.reset()
		samp.reset()

	button.on_clicked(reset)

	rax = plt.axes([0.025, 0.5, 0.15, 0.15], facecolor=axcolor)
	radio = RadioButtons(rax, ('red', 'blue', 'green'), active=0)

	def colorfunc(label):
		l.set_color(label)
		fig.canvas.draw_idle()

	radio.on_clicked(colorfunc)

	plt.show()
def merge_DSRC(experiment_path,filename,labels=None):
	###################### DSRC####################
	count = 0

	df_merged = []
	for day in os.listdir(experiment_path):
		if not 'March' in day:
			break
		for experiment_num_1 in os.listdir(os.path.join(experiment_path, day)):
			for bs in os.listdir(os.path.join(experiment_path, day, experiment_num_1)):
				for experiment_num_2 in os.listdir(os.path.join(experiment_path, day, experiment_num_1,bs)):
					path = os.path.join(experiment_path, day, experiment_num_1,bs, experiment_num_2)

					df_list = DSCR.get_DSRC_dataframes(path) if labels is None else DSCR.get_DSRC_dataframes(path,labels)
					# print(path)
					result_list = []
					for idx, df in enumerate(df_list):
						if count == -1:
							break
						count += 1
						df = pd.DataFrame(df)
						# append Vehicle ID
						df['cpid'] = df['UniqueOBE_ID_Alias'].apply(lambda x: int(int(re.sub("[^0-9]", "", x)) - 51))
						result_list += [df]
					df_merged = [pd.concat(df_merged+result_list,axis=0)]
				for x in result_list:del x
				del result_list

			# append distance and angle only for RXE
	df = df_merged[0].sort_values('TimeStamp_ms', ascending=True)
	df.to_csv(os.path.join(experiment_path, '{}.csv'.format(filename)), index=False)
	return df


def add_seq_mark(path):
	merge_DSRC(path,['TimeStamp_ms','MsgSeqNum'])
def nonlinear_cluster(data,k):
	cluster = cl.AgglomerativeClustering(n_clusters=k)
	res = cluster.fit(np.array(data).reshape((-1, 1)))
	# cluster.cluster_centers_ = np.sort(cluster.cluster_centers_, axis=0)
	return res
def kmeans(data,k):
	cluster = cl.KMeans(n_clusters=k, init='k-means++', n_init=10, max_iter=300, tol=0.0001, precompute_distances='auto', verbose=0,
	                                 random_state=None, copy_x=True, n_jobs=1)
	res = cluster.fit(np.array(data).reshape((-1,1)))
	cluster.cluster_centers_=np.sort(cluster.cluster_centers_,axis=0)
	return res
def dummy(samples):
	x = samples['Range'].values
	y = samples['angle'].values

	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.scatter(samples['Range'].values, samples['angle'].values)
	ax.set_xlabel('Range')
	ax.set_ylabel('Distance')
	plt.show()
def heatmap_polar(title,x,y,z,proj=None):
	fig = plt.figure()
	# ax = fig.add_subplot(211)
	ax2 = fig.add_subplot(111, projection='polar')
	plt.title(title)
	plt.scatter(y,x,c=z,cmap='plasma',s=.05)
	plt.colorbar()
	plt.show(block=False)
	return



def get_same_direction_dataset(samples):
	samples_1 = samples[samples.ego_lane < 8]
	samples_2 = samples[samples.ego_lane > 7]
	samples_1_1 = samples_1[samples_1.contact_lane < 8]
	samples_2_2 = samples_2[samples_2.contact_lane > 7]
	res = pd.concat([samples_1_1,samples_2_2],axis=0)
	return res


def get_opposite_direction_dataset(samples):
	samples_1 = samples[samples.ego_lane < 8]
	samples_2 = samples[samples.ego_lane > 7]
	samples_1_2 = samples_1[samples_1.contact_lane> 7]
	samples_2_1 = samples_2[samples_2.contact_lane < 8]
	res = pd.concat([samples_1_2, samples_2_1], axis=0)
	return res
def plot_3d(title, x, y, z,xtitle,ytitle,ztitle,cmap='plasma'):
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	ax.scatter(x, y,z,c=z, cmap=cmap,s=.01)
	ax.set_title(title)
	plt.show(block=False)



def plot_surface(title, x, y, z):
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	ax.scatter(x, y, z, c=z, cmap='plasma', s=.01)
	ax.set_title(title)
	plt.show(block=False)
def rotate_df(df, theta,  lat_center, long_center,in_labels = ['long','lat']):
	T = np.array([[long_center], [lat_center]])
	P_T = np.array([df[in_labels[0]].values, df[in_labels[1]].values])-T
	rotation_matrix = np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])
	ll_r = np.matmul(rotation_matrix,P_T)+T
	df[in_labels[0]+'_rec']=ll_r[0,:]
	df[in_labels[1]+'_rec'] = ll_r[1, :]
	return df
def rotate(lat_long,theta,title_long,title_lat,lat_center,long_center):
	res = {}
	y = lat_long[0]-lat_center
	x = lat_long[1]-long_center
	res[title_long] = ((np.cos(theta)*x)-(np.sin(theta)*y))+long_center
	res[title_lat] = ((np.sin(theta) * x) + (np.cos(theta)* y))+lat_center
	# res[titles_long] = lat_long[1]
	# res[title_lat] = lat_long[0]
	return pd.Series(res)


def rectify_latlong_curved_fast(df, param_dict=None):
	if param_dict is None:
		df_param = df.sample(10000)
		df_start_lat = df_param[df_param.long < -118.069].lat.min()
		df_end_lat = df_param[df_param.long > -118.053].lat.min()
		df_start_long = df_param.long.min()
		long_dif = df_param.long.max() - df_start_long
		lat_dif = df_end_lat - df_start_lat
		theta = -np.arctan(lat_dif / long_dif)
	else:
		theta = param_dict['theta']
		df_start_lat = param_dict['lat_center']
		df_start_long = param_dict['long_center']
	df = rotate_df(df, theta, df_start_lat, df_start_long)
	df = rotate_df(df, theta, df_start_lat, df_start_long, in_labels=['rxlong', 'rxlat'])

	return df
def rectify_latlong_fast(df, param_dict=None):
	if param_dict is None:
		df_param = df.sample(1000)
		df_start_lat = df_param[df_param.long < -118.069].lat.min()
		df_end_lat = df_param[df_param.long > -118.053].lat.min()
		df_start_long = df_param.long.min()
		long_dif = df_param.long.max() - df_start_long
		lat_dif = df_end_lat - df_start_lat
		theta = -np.arctan(lat_dif / long_dif)
	else:
		theta = param_dict['theta']
		df_start_lat = param_dict['lat_center']
		df_start_long = param_dict['long_center']
	df = rotate_df(df,theta,df_start_lat,df_start_long)
	df = rotate_df(df, theta, df_start_lat, df_start_long,in_labels=['rxlong','rxlat'])
	return df
def rectify_latlong(df,param_dict=None):
	if param_dict is None:
		df_start_lat = df[df.long < -118.069].lat.min()
		df_end_lat = df[df.long>-118.053].lat.min()
		df_start_long = df.long.min()
		long_dif = df.long.max() - df_start_long
		lat_dif = df_end_lat - df_start_lat
		theta = -np.arctan(lat_dif / long_dif)
	else:
		theta = param_dict['theta']
		df_start_lat= param_dict['lat_center']
		df_end_lat= param_dict['long_center']

	a = df[['lat', 'long']].apply(lambda x: rotate(x, theta, 'long_rec', 'lat_rec',lat_center=df_start_lat, long_center=df.long.min()), axis=1)
	b = df[['rxlat', 'rxlong']].apply(lambda x: rotate(x, theta, 'rxlong_rec', 'rxlat_rec', lat_center=df_end_lat, long_center=df.long.min()),
	                                  axis=1)
	df = pd.concat([df, a, b], axis=1)
	param_dict = {'lat_center':df_start_lat,'long_center':df_start_long,'theta':theta}
	return df
def crop_range(df):
	df = df[df.long >= DSCR.SEALBEACH_COORD['long_min']]
	df = df[df.long < DSCR.SEALBEACH_COORD['long_max']]
	df = df[df.rxlong >= DSCR.SEALBEACH_COORD['long_min']]
	df = df[df.rxlong < DSCR.SEALBEACH_COORD['long_max']]
	return df
def set_lane(samples):
	cluster = kmeans(samples['rxlat_rec'], 16)
	lat = np.array(samples['lat_rec'].values).reshape(-1, 1)
	rxlat = np.array(samples['rxlat_rec'].values).reshape(-1, 1)
	samples['contact_lane'] = cluster.predict(lat)
	samples['ego_lane'] = cluster.predict(rxlat)
	return samples


def set_lane_curved_contact(samples):
	# cluster = nonlinear_cluster(samples[['rxlat_rec','rxlong_rec']], 2)
	# subsamples = samples.sample(n=100000)
	# cluster = cl.Birch(n_clusters=16)
	cluster = cl.KMeans(n_clusters=16, n_init=20, tol=10)
	# cluster = cl.DBSCAN()
	# cluster = cl.Birch(n_clusters=16)
	total_batch = 10
	max_long = samples.long_rec.max()
	min_long = samples.long_rec.min()
	step = (max_long - min_long) / total_batch
	res = []
	for batch_index in np.arange(total_batch):
		start_val = (batch_index * step) + min_long
		end_val = ((batch_index + 1) * step) + min_long
		print(start_val, end_val)
		temp_samples = samples[samples.long_rec > start_val]
		temp_samples = temp_samples[temp_samples.long_rec < end_val]
		lat = np.transpose(np.array([temp_samples['lat_rec'].values * 1000000000, temp_samples['long_rec'].values / 100]))
		# lat_fit = np.transpose(np.array([subsamples['lat_rec'].values * 1000, subsamples['long_rec'].values]))
		# 	cluster =cluster.partial_fit(lat)
		# 	cluster=cluster.partial_fit(rxlat)
		cluster.fit(lat)
		centers = cluster.cluster_centers_
		cluster.cluster_centers_ = centers[centers.argsort(0)[:, 0]]
		temp_samples['contact_lane'] = cluster.predict(lat)
		res += [temp_samples]
		print(batch_index , '/' , total_batch)
	res_concat = pd.concat(res)
	res_concat = res_concat.sort_values('TimeStamp_ms', ascending=True)
	return res_concat
def set_lane_curved_ego(samples):
	# cluster = nonlinear_cluster(samples[['rxlat_rec','rxlong_rec']], 2)
	# subsamples = samples.sample(n=100000)
	# cluster = cl.Birch(n_clusters=16)
	cluster = cl.KMeans(n_clusters=16,n_init=10,tol=1)
	# cluster = cl.DBSCAN()
	# cluster = cl.Birch(n_clusters=16)
	total_batch = 10
	max_long = samples.long_rec.max()
	min_long = samples.long_rec.min()
	step = (max_long-min_long)/total_batch
	res = []
	for batch_index in np.arange(total_batch):
		start_val = (batch_index*step)+min_long
		end_val = ((batch_index+1) * step) + min_long
		print(start_val,end_val)
		temp_samples = samples[samples.rxlong_rec > start_val]
		temp_samples = temp_samples[temp_samples.rxlong_rec <end_val]
		rxlat = np.transpose(np.array([temp_samples['rxlat_rec'].values*1000000000,temp_samples['rxlong_rec'].values/100]))
	# lat_fit = np.transpose(np.array([subsamples['lat_rec'].values * 1000, subsamples['long_rec'].values]))
	# 	cluster =cluster.partial_fit(lat)
	# 	cluster=cluster.partial_fit(rxlat)
		cluster.fit(rxlat)
		centers = cluster.cluster_centers_
		cluster.cluster_centers_ = centers[centers.argsort(0)[:,0]]
		temp_samples['ego_lane']=cluster.predict(rxlat)
		res+=[temp_samples]
		print(batch_index,'/',total_batch)
		# print(len(cluster.cluster_centers_))
	res_concat = pd.concat(res)
	res_concat = res_concat.sort_values('TimeStamp_ms', ascending=True)
	return res_concat
def concat_distance_angle(df):

	# df= df.drop(['distance','angle'],axis=1)
	df = pd.concat([df, df[['rxlat_rec', 'rxlong_rec', 'lat_rec', 'long_rec']].apply(DSCR.get_polar_list_input, axis=1)], axis=1)
	a = df[df.ego_lane > 7]['angle'].values
	df.loc[df.ego_lane > 7,'angle'] =  - a
	return df

def get_density(timestampms, den):
	bias = 1490511600
	timestamp_unbiased=(timestampms//1000 - bias)
	timestamp = ((timestamp_unbiased//300)*300)+bias
	# print(timestamp)
	# time = den[den.TimeStamp_S < timestamp[0]].TimeStamp_S.max()
	density = den[den.TimeStamp_S == timestamp[0]]
	density = density[['1', '2', '3', '4', '5', '6', '7', '8', 'Total']]
	density=density.to_dict(orient='list')
	return pd.Series(density)
def add_timestamp_s_5min(df):
	bias = 1490511600
	df['TimeStamp_S']=((((df.TimeStamp_ms//1000)- bias)//300)*300)+bias
	return df
def concat_density_optimized(df):

	we_den_df = pd.read_csv('../Dataset/DSRC/west.csv')
	ea_den_df = pd.read_csv('../Dataset/DSRC/east.csv')
	density_df =pd.DataFrame()
	ea_den_df.rename(columns={'Total':'Total_East_Density','1':'7','2':'6','3':'5','4':'4','5':'3','6':'2','7':'1','8':'0'},inplace=True)
	we_den_df.rename(columns={'Total': 'Total_West_Density', '1': '8', '2': '9', '3': '10', '4': '11', '5': '12', '6': '13', '7': '14', '8': '15'},
	                 inplace=True)
	den_df = pd.merge(ea_den_df,we_den_df,on='TimeStamp_S')
	den_df['Av_East_Density'] = den_df['Total_East_Density']/8
	den_df['Av_West_Density'] = den_df['Total_West_Density'] / 8
	df = add_timestamp_s_5min(df)
	q = pd.merge(df,den_df,on='TimeStamp_S')
	return q
def parallel_df(df,func):
	df_split = np.array_split(df, 10)
	pool = Pool(8)
	df = pd.concat(pool.map(func, df_split))
	pool.close()
	pool.join()
	return df

def concat_density(df):
	def multiply_columns(df,den_df):
		den_merged_df = df[['TimeStamp_ms']].apply(lambda x: get_density(x, den_df), axis=1)
		return den_merged_df
	we_den_df = pd.read_csv('../Dataset/DSRC/west.csv')
	ea_den_df = pd.read_csv('../Dataset/DSRC/east.csv')
	# den_east = parallel_df(den_east,)
	den_east = df[['TimeStamp_ms']].apply(lambda x:get_density(x,ea_den_df),axis=1)
	print('east done')
	den_west = df[['TimeStamp_ms']].apply(lambda x: get_density(x, we_den_df),axis=1)
	print('west done')
	den_west.columns=['8', '9', '10', '11', '12','13','14','15','total_west_density']
	den_east.columns =['7','6','5','4','3','2','1','0','total_east_density']
	# den_west_global = den_west.rename({'1':'8', '2': '9', '3': '10', '4': '11', '5': '12', '6': '13', '7': '14','8':'15',
	#                                    'Total':'total_west_density'})
	# den_east_global = den_east.rename({'1': '7', '2': '6', '3': '5', '4': '4', '5': '3', '6': '2', '7': '1', '8': '0', 'Total': 'total_east_density'})
	den = pd.concat([den_west,den_east],axis=1)
	den['Total_density'] = den[['total_east_density','total_west_density']].apply(lambda x: x[0][0]+x[1][0],axis=1)
	df = pd.concat([df,den],axis=1)
	return df

def prepare_dataset(df):
	df = crop_range(df)
	print('crop done')
	df = rectify_latlong_fast(df)
	print('rotate done')
	df = set_lane_curved_ego(df)
	print('ego_done')
	df = set_lane_curved_contact(df)
	print('lane id done')
	df = concat_distance_angle(df)
	print('polar done')
	df = concat_density_optimized(df)
	print('density done')
	return df
def plot_density_distance_angle_RSS_lanes_with_density_sep(df):
	tempdf = df
	density_range = [tempdf.Total_density.min(), 400]
	density_length = density_range[1]-density_range[0]
	split = 3
	for e_lane in np.arange(16):
		tempdf_ego=df[df.ego_lane==e_lane]
		tempdf_ego_c = tempdf_ego
		for idx,density in enumerate(np.linspace(density_range[0],density_range[1],split)):
			tempdf = tempdf_ego_c[tempdf_ego_c.Total_density>=density]
			r2 = density + (density_length / split)
			tempdf = tempdf[tempdf.Total_density < r2]
			tempdf = tempdf[tempdf.distance < 400]
			angle = tempdf['angle'].values;
			y = tempdf['distance'].values;
			x = tempdf['RSS'].values

			heatmap('ZZZ-RSS-Range-Angle_Distance_lane={}_density{}-{}'.format(e_lane,int(density),int(r2)), y, angle, x,
			        'plasma', '',
			        'distance', 'RSS')


def running_stat(X, Y, total_bins, fun_handle):
	bins = np.linspace(X.min(), X.max(), total_bins)
	delta = bins[1] - bins[0]
	centers = []
	idx = np.digitize(X, bins)
	running_median = []
	for k in np.arange(total_bins):
		val_points = Y[idx == k]
		if len(val_points) == 0:
			continue
		else:
			running_median += [fun_handle(Y[idx == k])]
			centers += [bins[k]]
	return centers - delta / 2, running_median


def running_stat3d(X1,X2, Y, total_bins, fun_handle):
	binsX1 = np.linspace(X1.min(), X1.max(), total_bins)
	deltaX1 = binsX1[1] - binsX1[0]
	centersX1 = []
	idxX1 = np.digitize(X1, binsX1)
	binsX2 = np.linspace(X2.min(), X2.max(), total_bins)
	deltaX2 = binsX2[1] - binsX2[0]
	centersX2 = []
	idxX2 = np.digitize(X2, binsX2)
	running_median = []
	for x1 in np.arange(total_bins):
		for x2 in np.arange(total_bins):
			val_points = Y[(idxX1 == x1) and(idxX2==x2)]
			if len(val_points) == 0:
				continue
			else:
				running_median += [fun_handle(val_points)]
				centersX1 += [binsX1[x1]]
				centersX2 += [binsX2[x2]]
	return centersX1 - deltaX1 / 2, running_median, centersX2 - deltaX2 / 2
def plot_different_density_RSS_distance(df, dir, show):
	def agg_plot(df_list, plot_group_name, dir):
		if not os.path.exists(dir):
			os.mkdir(dir)
		for i in ['Min', 'Max', 'Mean']:
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
			ego += [df['ego_lane'].mean()]
		# Polar all_range
		multiline_plot('{} Density AGG RSS-Distance Mean'.format(plot_group_name), x_range_list=range, y_theta_list=RSS, show=show,
		               dir=path.join(dir, 'Mean'), xtitle='Range(meter)', ytitle='RSS', fun_handle=np.mean, line_id_list=ego)
		multiline_plot('{} Density AGG RSS-Distance Median'.format(plot_group_name), x_range_list=range, y_theta_list=RSS, show=show,
		               dir=path.join(dir, 'Median'), xtitle='Range(meter)', ytitle='RSS', fun_handle=np.median, line_id_list=ego)
		multiline_plot('{} Density AGG RSS-Distance Max'.format(plot_group_name), x_range_list=range, y_theta_list=RSS, show=show,
		               dir=path.join(dir, 'Max'), xtitle='Range(meter)', ytitle='RSS', fun_handle=np.max, line_id_list=ego)
		multiline_plot('{} Density AGG RSS-Distance Min'.format(plot_group_name), x_range_list=range, y_theta_list=RSS, show=show,
		               dir=path.join(dir, 'Min'), xtitle='Range(meter)', ytitle='RSS', fun_handle=np.min, line_id_list=ego)

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
	for dif in np.arange(16):
		diff_df = []
		for e_lane in np.arange(16):
			lane_dir = os.path.join(dir, str(e_lane))
			tempdf_ego = df[df.ego_lane == e_lane]
			tempdf_ego = tempdf_ego[np.abs(tempdf_ego.ego_lane - tempdf_ego.contact_lane) == dif]
			diff_df += [tempdf_ego]
			plotter(tempdf_ego, 'Lane#{} Lane Difference:{}'.format(e_lane, dif), os.path.join(dir, 'Lane Diff:' + str(dif)))
		agg_plot(diff_df, 'Lane Difference:{}'.format(dif), os.path.join(dir))


def plot__RSS_distance_lane_group_multiplot_density_seperate(df, dir, show, lane_set_cardinality=1, mirror_lane_numbering=False,density_bins=5):
	def agg_plot(df_list, plot_group_name, dir):
		if not os.path.exists(dir):
			os.mkdir(dir)
		for i in ['Min', 'Max', 'Mean', 'Median','Samples']:
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
		multiline_plot('{}RSS-Distance Samples'.format(plot_group_name), x_range_list=range, y_theta_list=RSS, show=show, dir=path.join(dir,
		                                                                                                                                'Samples'),
		               xtitle='Range(meter)', ytitle='Samples', fun_handle=len, line_id_list=ego)

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
	min_density  = df.Total_density.min()
	step_density = (max_density-min_density)/density_bins
	density_list = min_density+(step_density*np.arange(density_bins+1))
	for i in np.arange(density_bins):
		diff_df = []
		for e_lane in np.arange(16 / lane_set_cardinality):
			lane_dir = os.path.join(dir, str(e_lane))
			if mirror_lane_numbering:
				df_east = df[df.ego_lane < 8]
				df_west = df[df.ego_lane >= 8]

				tempdf_ego = pd.concat(
					[df_east[(df_east.ego_lane // lane_set_cardinality) == e_lane], df_west[((15 - df_west.ego_lane) // lane_set_cardinality) ==
					                                                                        e_lane]])
			else:
				tempdf_ego = df[(df.ego_lane // lane_set_cardinality) == e_lane]
			tempdf_ego = tempdf_ego[tempdf_ego.Total_density<density_list[i+1]]
			tempdf_ego = tempdf_ego[tempdf_ego.Total_density > density_list[i]]
			diff_df += [tempdf_ego]
		# plotter(tempdf_ego, 'Lane#{} Lane Difference:{}'.format(e_lane, dif), os.path.join(dir, 'Lane Diff:' + str(dif)))
		agg_plot(diff_df, 'DensityRange {}:{} '.format(int(density_list[i]/16),int(density_list[i+1]/16)), os.path.join(dir))
def plot__RSS_distance_lane_group_multiplot(df, dir, show, lane_set_cardinality=1, mirror_lane_numbering=False):
	def agg_plot(df_list, plot_group_name, dir):
		if not os.path.exists(dir):
			os.mkdir(dir)
		for i in ['Min', 'Max', 'Mean', 'Median','Samples']:
			if not path.exists(path.join(dir, i)):
				os.mkdir(path.join(dir, i))
		theta = []
		range = []
		RSS = []
		density = []
		ego = []
		for df in df_list:
			range += [df['distance'].values]
			RSS += [df['RSS'].values]
			density += [df['Total_density'].values]
			ego += [set(df['ego_lane'].values)]
		# Polar all_range
		multiline_plot('{}RSS-Distance Mean'.format(plot_group_name), x_range_list=range, y_theta_list=RSS, show=show,
		               dir=path.join(dir, 'Mean'), xtitle='Range(meter)', ytitle='RSS', fun_handle=np.mean, line_id_list=ego)
		multiline_plot('{}RSS-Distance Median'.format(plot_group_name), x_range_list=range, y_theta_list=RSS, show=show,
		               dir=path.join(dir, 'Median'), xtitle='Range(meter)', ytitle='RSS', fun_handle=np.median, line_id_list=ego)
		multiline_plot('{}RSS-Distance Max'.format(plot_group_name), x_range_list=range, y_theta_list=RSS, show=show,
		               dir=path.join(dir, 'Max'), xtitle='Range(meter)', ytitle='RSS', fun_handle=np.max, line_id_list=ego)
		multiline_plot('{}RSS-Distance Min'.format(plot_group_name), x_range_list=range, y_theta_list=RSS, show=show,
		               dir=path.join(dir, 'Min'), xtitle='Range(meter)', ytitle='RSS', fun_handle=np.min, line_id_list=ego)
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
		# tempdf_ego = tempdf_ego[np.abs(tempdf_ego.ego_lane - tempdf_ego.contact_lane) == dif]
		diff_df += [tempdf_ego]
		# plotter(tempdf_ego, 'Lane#{} Lane Difference:{}'.format(e_lane, dif), os.path.join(dir, 'Lane Diff:' + str(dif)))
	agg_plot(diff_df, '', os.path.join(dir))
def plot_lane_difference_RSS_distance(df,dir,show,lane_set_cardinality=1,mirror_lane_numbering=False):
	def agg_plot(df_list, plot_group_name, dir):
		if not os.path.exists(dir):
			os.mkdir(dir)
		for i in ['Min','Max','Mean','Median','Samples']:
			if not path.exists(path.join(dir,i)):
				os.mkdir(path.join(dir,i))
		theta=[]
		range =[]
		RSS =[]
		density = []
		ego = []
		for df in df_list:
			theta += [df['angle'].values]
			range += [df['distance'].values]
			RSS += [df['RSS'].values]
			density += [df['Total_density'].values]
			ego += [set(df['ego_lane'].values)]
		# Polar all_range
		multiline_plot('{} RSS-Distance Mean'.format(plot_group_name), x_range_list=range, y_theta_list=RSS,
		               show=show, dir=path.join(dir,'Mean'), xtitle='Range(meter)', ytitle='RSS',fun_handle=np.mean, line_id_list=ego)
		multiline_plot('{} RSS-Distance Median'.format(plot_group_name), x_range_list=range, y_theta_list=RSS, show=show, dir=path.join(
			dir, 'Median'),
		               xtitle='Range(meter)', ytitle='RSS', fun_handle=np.median, line_id_list=ego)
		multiline_plot('{} RSS-Distance Max'.format(plot_group_name), x_range_list=range, y_theta_list=RSS, show=show, dir=path.join(dir,
		                                                                                                                                         'Max'),
		               xtitle='Range(meter)', ytitle='RSS', fun_handle=np.max, line_id_list=ego)
		multiline_plot('{} RSS-Distance Min'.format(plot_group_name), x_range_list=range, y_theta_list=RSS, show=show, dir=path.join(dir,
		                                                                                                                                         'Min'),
		               xtitle='Range(meter)', ytitle='RSS', fun_handle=np.min, line_id_list=ego)
		multiline_plot('{} RSS-Distance Samples'.format(plot_group_name), x_range_list=range, y_theta_list=RSS, show=show,
		               dir=path.join(dir, 'Samples'), xtitle='Range(meter)', ytitle='Samples', fun_handle=len, line_id_list=ego)
		# multiline_plot('{} Density AGG RSS-Distance std'.format(plot_group_name), x_range_list=range, y_theta_list=RSS, show=show, dir=path.join(dir,
		#                                                                                                                                          'Mean'),
		#                xtitle='Range(meter)', ytitle='RSS', fun_handle=np.std, line_id_list=ego)
	def plotter(df, plot_group_name, dir):
		if not os.path.exists(dir):
			os.mkdir(dir)
		theta = df['angle'].values
		range = df['distance'].values
		RSS = df['RSS'].values
		density = df['Total_density'].values
		# Polar all_range
		heatmap('{} Density Heat Map RSS-Distance'.format(plot_group_name), x_range=range, y_theta=RSS, color_param=density,
		        proj='',
		        show=show, dir=dir, xtitle='Range(meter)', ytitle='RSS', colorbar_label='density(dbm)')
	if not os.path.exists(dir):
		os.mkdir(dir)
	for dif in np.arange(16):
		diff_df = []
		for e_lane in np.arange(16/lane_set_cardinality):
			lane_dir = os.path.join(dir, str(e_lane))
			if not os.path.exists(lane_dir):
				os.mkdir(lane_dir)
			if mirror_lane_numbering:
				df_east = df[df.ego_lane<8]
				df_west = df[df.ego_lane>=8]

				tempdf_ego = pd.concat([df_east[(df_east.ego_lane//lane_set_cardinality) == e_lane], df_west[
					((15-df_west.ego_lane) // lane_set_cardinality) == e_lane]])
			else:
				tempdf_ego = df[(df.ego_lane // lane_set_cardinality) == e_lane]
			tempdf_ego = tempdf_ego[np.abs(tempdf_ego.ego_lane-tempdf_ego.contact_lane)==dif]
			diff_df+=[tempdf_ego]
			plotter(tempdf_ego, 'Lane#{} Lane Difference:{}'.format(e_lane,dif), os.path.join(dir,'Lane Diff:'+str(dif)))
		agg_plot(diff_df, 'Lane Difference: {}'.format( dif), os.path.join(dir))
		# plotter(pd.concat((tempdf),axis=0),'Lane Difference:{}'.format(dif),dir)
def plot_density_distance_angle_RSS_lanes(df,dir,show):
	def plotter(df, plot_group_name, dir):
		if not os.path.exists(dir):
			os.mkdir(dir)
		theta = df['angle'].values
		range = df['distance'].values
		RSS = df['RSS'].values
		density = df['Total_density'].values
		# Polar all_range
		heatmap('{} RSS Heat Map Range-Theta in Cartisian'.format(plot_group_name), x_range=range, y_theta=theta, color_param=RSS, proj='polar',
		        show=show, dir=dir, xtitle='Range(meter)', ytitle='Theta', colorbar_label='RSS(dbm)')
		heatmap('{} RSS Heat Map Range-Theta in Polar'.format(plot_group_name), x_range=range, y_theta=theta, color_param=RSS, proj='', show=show,
		        dir=dir, xtitle='Range(meter)', ytitle='Theta', colorbar_label='RSS(dbm)')
		heatmap('{} Density Heat Map Range-RSS'.format(plot_group_name), x_range=range, y_theta=RSS, color_param=density, proj='', show=show,
		        dir=dir,
		        xtitle='Range(meter)', ytitle='RSS(dbm)', colorbar_label='Density')
		heatmap('{} RSS Heat Map Range-Density'.format(plot_group_name), x_range=density, y_theta=range, color_param=RSS, proj='', show=show,
		        dir=dir,
		        ytitle='Range(meter)', xtitle='Density', colorbar_label='RSS(dbm)')
		plot3d('{} RSS Range-Theta 3D in Polar'.format(plot_group_name), x_range=range, y_theta=theta, value=RSS, color_param=RSS, proj='',
		       show=show,
		       dir=dir, xtitle='Range', ytitle='Theta', ztitle='RSS(dbm)', colorbar_label='RSS(dbm)')
		plot3d('{} RSS Range-Theta 3D in Cartesiean'.format(plot_group_name), x_range=np.cos(theta) * range, y_theta=np.sin(theta) * range,
		       value=RSS,
		       color_param=RSS, show=show, dir=dir, xtitle='Road Direction', ytitle='Road Width Direction', ztitle='RSS(dbm)',
		       colorbar_label='RSS(dbm)')
		print(plot_group_name)
	if not os.path.exists(dir):
		os.mkdir(dir)
	for e_lane in np.arange(16):
		lane_dir = os.path.join(dir,str(e_lane))
		if not os.path.exists(lane_dir):
			os.mkdir(lane_dir)
		tempdf_ego = df[df.ego_lane == e_lane]
		plotter(tempdf_ego,'Lane{}'.format(e_lane),lane_dir)
def realdomain_to_descreet(df,field_list,n):
	df_res = pd.DataFrame()
	df = df.loc[:, ~df.columns.str.contains('^Unnamed')]
	for field in df.columns:
		print(field)
		if field in field_list:
			max_val = df[field].max()
			min_val = df[field].min()
			range = max_val-min_val
			bin_range = range/n
			df_res[field] = (((df[field]-min_val)//bin_range)*bin_range) +min_val
		else:
			df_res[field]=df[field]
	return df_res

def plot_distance_angle_RSS(df,show=True,dir=''):
	def get_range_theta_RSS(df,max_range=-1):
		res = df if max_range==-1 else df[df.distance<max_range]
		return res['distance'].values,res['angle'].values,res['RSS'].values,res['Total_density'].values
	# # RSS VS angle VS distance
	def plotter(df,plot_group_name,dir):
		if not os.path.exists(dir):
			os.mkdir(dir)
		theta= df['angle'].values
		range = df['distance'].values
		RSS = df['RSS'].values
		density = df['Total_density'].values
		# Polar all_range
		heatmap('{} RSSI Heat Map Range-Theta in Cartisian'.format(plot_group_name),x_range=range,y_theta=theta,color_param=RSS,proj='polar',
		        show=show,dir=dir,
		        xtitle='Range(meter)',ytitle='Theta',colorbar_label='RSS(dbm)')
		heatmap('{} RSSI Heat Map Range-Theta in Polar'.format(plot_group_name), x_range=range, y_theta=theta, color_param=RSS, proj='', show=show,
		        dir=dir,
		        xtitle='Range(meter)', ytitle='Theta', colorbar_label='RSS(dbm)')
		heatmap('{} Density Heat Map Range-RSSI'.format(plot_group_name), x_range=range, y_theta=RSS, color_param=density, proj='', show=show,
		        dir=dir,
		        xtitle='Range(meter)', ytitle='RSS(dbm)', colorbar_label='Density')
		heatmap('{} RSSI Heat Map Range-Density'.format(plot_group_name), x_range=density, y_theta=range, color_param=RSS, proj='', show=show,
		        dir=dir,
		        ytitle='Range(meter)', xtitle='Density', colorbar_label='RSS(dbm)')
		plot3d('{} RSSI Range-Theta 3D in Polar'.format(plot_group_name), x_range=range, y_theta=theta,value=RSS, color_param=RSS, proj='', show=show,
		       dir=dir,
		       xtitle
		='Range',
		       ytitle='Theta',ztitle='RSS(dbm)', colorbar_label = 'RSS(dbm)')
		plot3d('{} RSSI Range-Theta 3D in Cartesiean'.format(plot_group_name), x_range=np.cos(theta)*range, y_theta=np.sin(theta)*range, value=RSS,
		       color_param=RSS,
		       show=show,
		       dir=dir, xtitle='Road Direction', ytitle='Road Width Direction', ztitle='RSS(dbm)', colorbar_label='RSS(dbm)')
		print(plot_group_name)
	plotter(df,'All Lanes',os.path.join(dir,'both direction'))
	plotter(get_same_direction_dataset(df), 'Same Direction', os.path.join(dir, 'same direction'))
	plotter(get_opposite_direction_dataset(df), 'Opposite Direction', os.path.join(dir, 'oposite direction'))
def multiline_plot(title, x_range_list, y_theta_list, line_id_list, xtitle='', ytitle='', dir='', show=True, save=True,
                   fun_handle=np.mean):
	fig = plt.figure(figsize=(15, 10))
	ax2 = fig.add_subplot(111)
	plt.title(title)
	count=0
	for x,y,id in zip(x_range_list,y_theta_list,line_id_list):
		count+=1
		if len(x)==0 or len(y)==0:
			print(x)
			continue
		x_bin,y_med= running_stat(x,y,75,fun_handle)
		plt.plot(x_bin, y_med,'--',label = 'Lane:{}'.format(id))

		print(count)
	# plt.hold(True)
	plt.legend()
	plt.xlabel(xtitle)
	plt.ylabel(ytitle)
	if (not dir == '') and save:
		plt.savefig(os.path.join(dir, title + '.png'))
	plt.show(block=False) if show else plt.close()

def heatmap(title, x_range, y_theta, color_param, cmap='plasma', proj=None, xtitle='', ytitle='',dir='',show=True,save=True,s=.1,colorbar_label=''):
	fig = plt.figure(figsize=(15,10))
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


def plot3d(title, x_range, y_theta,value, color_param, cmap='plasma', proj=None, xtitle='', ytitle='',ztitle='', dir='', show=True, save=True, s=.1,
           colorbar_label=''):
	fig = plt.figure(figsize=(15, 10))
	ax = fig.add_subplot(111, projection='3d')
	p = ax.scatter(x_range, y_theta, value, '.', c=color_param, cmap=cmap, s=.05)
	ax.set_xlabel(xtitle)
	ax.set_ylabel(ytitle)
	ax.set_zlabel(ztitle)
	cb = fig.colorbar(p)
	cb.set_label(colorbar_label)
	plt.title(title)
	for idx,angle in enumerate(np.linspace(0, 90,3)):
		ax.view_init(0, angle)
		plt.draw()
		plt.pause(.001)
		if (not dir == '') and save:
			plt.savefig(os.path.join(dir, title +str(idx)+ 'view0.png'))
		# plt.show(block=False) if show else plt.close()
	for idx, angle in enumerate(np.linspace(0, 90, 3)):
		ax.view_init(30, angle)
		plt.draw()
		plt.pause(.001)
		if (not dir == '') and save:
			plt.savefig(os.path.join(dir, title + str(idx) + 'view30.png'))
		# plt.show(block=False) if show else plt.close()
	plt.show(block=False) if show else plt.close()
def fix_density(df):
	for i in ['total_east_density','total_west_density']:
		df[str(i)]=df[str(i)].apply(lambda x : float(x.replace('[','').replace(']','')))
	return df
def save_based_on_lane(df):
	df = crop_range(df)
	print('crop done')
	df = rectify_latlong_fast(df)
	print('rotate done')
def get_params(df):
	df = crop_range(df)
	print('crop done')
	df = rectify_latlong(df)
	print('rotate done')
	df = set_lane(df)
	print('lane id done')
	df = concat_distance_angle(df)
	print('polar done')
	df = concat_density(df)

if __name__ == '__main__':
	filename = 'all_merged.csv'
	experiment_path = '../Dataset/DSRC/'
	df_filepath = path.join(experiment_path,filename)
	plot_path = '../Plots/'
	field_list = [ 'lat', 'rxlat', 'rxlong', 'Range', 'RSS', 'long',

	               'lat_rec', 'long_rec', 'rxlat_rec', 'rxlong_rec', 'ego_lane', 'angle', 'distance', '8', '9', '10', '11', '12',
	               '13',
	               '14', '15', 'total_west_density', '7', '6', '5', '4', '3', '2', '1', '0', 'total_east_density', 'Total_density']
	df= pd.read_csv(df_filepath,index_col=False)
	# df = df.sample(n=4000)
	df = prepare_dataset(df)
	df.to_csv(path.join(experiment_path,'alldone_density.csv'))
	# heatmap('', x_range=df['rxlong_rec'].values, y_theta=df['rxlat_rec'].values, color_param=df['ego_lane'].values,  show=True,save=False,
     #    dir=dir, xtitle='Range(meter)', ytitle='Theta', colorbar_label='RSS(dbm)',cmap='tab20')
	heatmap('', x_range=df['long_rec'].values, y_theta=df['lat_rec'].values, color_param=df['contact_lane'].values, show=True, save=False, dir=dir,
	        xtitle='Range(meter)', ytitle='Theta', colorbar_label='RSS(dbm)', cmap='tab20')
	# plt.pause(1000)
# df_fixed = fix_density(df)
# 	df = realdomain_to_descreet(df,field_list,50)
# 	print(df)
# 	for cardinality in 2**np.arange(4):
# 		f1 = path.join(plot_path, 'RSS_Distance_lane_group_density_sep')
# 		f2 = path.join(plot_path, 'RSS_Distance_lane_group')
# 		if path.exists(f1)==-1: os.mkdir(f1)
# 		if path.exists(f2) == -1: os.mkdir(f2)
# 		plot__RSS_distance_lane_group_multiplot_density_seperate(df,path.join(plot_path,'RSS_Distance_lane_group_density_sep',
# 		                                                                      'Discrete_{}lane'.format(cardinality)),
# 		                                                         show=False,
# 		                                        lane_set_cardinality=cardinality,
# 		                                        mirror_lane_numbering=False,density_bins=5)
# 		plot__RSS_distance_lane_group_multiplot(df, path.join(plot_path, 'RSS_Distance_lane_group', 'Discrete_{}lane'.format(cardinality)),
# 	                                                         show=False, lane_set_cardinality=cardinality, mirror_lane_numbering=False)
#
# 		plot__RSS_distance_lane_group_multiplot_density_seperate(df, path.join(plot_path, 'RSS_Distance_lane_group_density_sep',
# 		                                                                       'Discrete_{}lane Mirror'.format(cardinality)), show=False,
# 		                                                         lane_set_cardinality=cardinality/2, mirror_lane_numbering=True, density_bins=5)
# 		plot__RSS_distance_lane_group_multiplot(df, path.join(plot_path, 'RSS_Distance_lane_group', 'Discrete_{}lane '
# 		                                                                                                        'Mirror'.format(cardinality)),
# 	                                        show=False, lane_set_cardinality=cardinality/2, mirror_lane_numbering=True)

	# merge_DSRC(experiment_path, 'merge1')

# 	plot_distance_angle_RSS(df,False, path.join(plot_path, 'RSS_angle_distance'))
# 	df.to_csv(path.join(experiment_path,'alldone1M_discrete50bin.csv'))
	# plot_lane_difference_RSS_distance(df_all, dir=os.path.join(plot_path, 'RSS_angle_distance_density'), show=False)
	# plt.pause(20)
	# for lane in np.arange(10):
	# 	prepare_dataset()
	# 	df = df_all[df_all.ego_lane==lane]
	# 	df.to_csv(os.path.join(experiment_path,'Lanes','lane{}.csv'))

	#
	# df = merge_DSRC(experiment_path)
	# df = pd.read_csv('../Dataset/DSRC/all_merged.csv')
	# df = df.sample(n=1000000)
	# plot_distance_angle_RSS(df)
	# plot_density_distance_angle_RSS(df)
	# sliding(df)
	# df = prepare_dataset(df)
	# df.to_csv(os.path.join(experiment_path,'alldone.csv'))
	# samples =df
	# samples = df[df.LogRecType == 'RXE'][['lat','rxlat','rxlong','LogRecType','Range', 'distance', 'angle','egoID','RSS','UniqueOBE_ID_Alias','long',
	#                                       'speed',
	#                                       'heading']]
	# plt.pause(20)
	# samples = samples.sample(n=10000)
	# samples.to_csv(os.path.join(experiment_path,'10k_merged.csv'))
	# del df
	# samples = samples[samples.distance<500 ]
	# samples = samples.sample(frac=.5)
	# x = samples['Range'].values;y = samples['angle'].values;z = samples['RSS'].values
	# cluster = kmeans(samples['rxlat_rec'],16)
	# lat = np.array(samples['lat_rec'].values).reshape(-1, 1)
	# rxlat = np.array(samples['rxlat_rec'].values).reshape(-1, 1)
	# samples['contact_lane'] = cluster.predict(lat)
	# samples['ego_lane'] = cluster.predict(rxlat)
	#
	# heatmap('lat_long',samples['lat_rec'],samples['long_rec'],samples['ego_lane'],cmap='Paired')
	# heatmap('reclat_long', samples['lat_rec'], samples['long_rec'], samples['ego_lane'], cmap='Paired')
	# heatmap('rxlat_long', samples['rxlat_rec'], samples['long_rec'], samples['lane'], cmap='Paired')
	# plot_3d('RSS-Range-Angle',x,y,z)
	# samples_temp = get_same_direction_dataset(samples)
	# x = samples_temp['Range'].values;
	# y = samples_temp['angle'].values;
	# z = samples_temp['RSS'].values
	# plot_3d('RSS-Range-Angle_same_dir', x, y, z)
	# samples_temp = get_opposite_direction_dataset(samples)
	# x = samples_temp['Range'].values;
	# y = samples_temp['angle'].values;
	# z = samples_temp['RSS'].values
	# plot_3d('RSS-Range-Angle_opposite_dir', x, y, z)
	# DONE

		# samples_temp = get_same_direction_dataset(samples)
		# x = samples_temp['Range'].values;
		# y = samples_temp['angle'].values;
		# z = samples_temp['RSS'].values
		# heatmap('RSS-Range-Angle_same_dir', x, y, z)
		# samples_temp = get_opposite_direction_dataset(samples)
		# x = samples_temp['Range'].values;
		# y = samples_temp['angle'].values;
		# z = samples_temp['RSS'].values
		# heatmap('RSS-Range-Angle_opposite_dir', x, y, z)


















