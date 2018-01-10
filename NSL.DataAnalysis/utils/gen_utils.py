import os
def experiment_decorator(fun,experiment_path):
	date = 'March 26'
	experiment_path = '../Dataset/DSRC'
	counter = 1
	# animate_gps(os.path.abspath('../Dataset/DSRC/{}/1/2/'.format(date)))
	for day in os.listdir(experiment_path):
		if not 'March' in day:
			break
		for experiment_num_1 in os.listdir(os.path.join(experiment_path, day)):
			for experiment_num_2 in os.listdir(os.path.join(experiment_path, day, experiment_num_1)):
				path = os.path.join(experiment_path, day, experiment_num_1, experiment_num_2)
				print(path)
				fun(path)