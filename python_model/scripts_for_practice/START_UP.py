#!/usr/bin/python
#START_UP.py

'''Get working environment set up to start working on the model.'''

import os

files = os.listdir('./scripts_for_practice/')
files = files[1:len(files)].append(files[0])
for f in os.listdir('./scripts_for_practice/'):
	if f.startswith('PRACTICE'):
		print(f)
    	dec =raw_input('\n\nExecute?')
    	if dec == 'y':
			execfile('./scripts_for_practice/%s' % f)
	else:
		pass

