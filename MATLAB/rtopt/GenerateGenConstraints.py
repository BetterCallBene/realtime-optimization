#!/usr/bin/python

import os
import tempfile

print 'Building - MATLAB - File'

pathTmp = tempfile.gettempdir() + '/%s'
files = [
		'tmpRTConstraintsFunction', 
		'tmpRTConstraintsJacobi', 
		'tmpRTConstraintsHesse', 
		'tmpRTConstraintsJacobiCount',
		'tmpRTConstraintsCount',
		'tmpRTConstraintsHesseCount',
		]

template = 'template_GenConstraints.m'
file_to = 'GenConstraints.m'

string_to_replace = [
					['q(1)', 'q(1, timestep)'],
					['q(2)', 'q(2, timestep)'],
					['q(3)', 'q(3, timestep)'],
					['q(4)', 'q(4, timestep)'],
					['v(1)', 'v(1, timestep)'],
					['v(2)', 'v(2, timestep)'],
					['v(3)', 'v(3, timestep)'],
					['omega(1)', 'omega(1, timestep)'],
					['omega(2)', 'omega(2, timestep)'],
					['omega(3)', 'omega(3, timestep)'],
					['u(1)', 'u(1, timestep)'],
					['u(2)', 'u(2, timestep)'],
					['u(3)', 'u(3, timestep)'],
					['u(4)', 'u(4, timestep)'],
					['*', '.*'],
					['/', './'],
					['^', u'.^'],
					['ones', 'onesV']
					]
print 'Temporaere loeschen'

for filename in files:
	path = pathTmp % filename

	if os.path.exists(path):
		os.remove(path)

if os.path.exists(file_to):
	os.remove(file_to)
	
print 'Do maple stuff'
os.system('maple < GenerateConstraints.mpl')
template_file = open(template, 'r')
file_write_to = open(file_to, 'w')

print 'Ersetzungen durchfuehren'
try:
	text = template_file.read()
finally:
	template_file.close()

i = 0

for filename in files:

	path = pathTmp % filename
	HFile = open(path, 'r')

	text_rep = u'$%d$' % i 

	
	try:
		text_input = HFile.read()
		for str_pair in string_to_replace:
			text_input = text_input.replace(str_pair[0], str_pair[1])
		
		text = text.replace(text_rep, text_input)
		
	finally:
		HFile.close()

	i = i + 1

file_write_to.write(text)
file_write_to.close()

print 'Temporaere Dateien loeschen'

for filename in files:
	path = pathTmp % filename

	if os.path.exists(path):
		os.remove(path)

 
