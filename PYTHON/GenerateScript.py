#!/usr/bin/python

import os
import tempfile

MAPLE_PATH = '/Library/Frameworks/Maple.framework/Versions/Current/bin/maple'

def generate(Output, TmpFiles, MplFile, TemplateFile, strings_to_replace):
	print 'Building - MATLAB - File'

	pathTmp = tempfile.gettempdir() + '/%s'

	print 'Temporaere loeschen'

	for filename in TmpFiles:
		path = pathTmp % filename[0]

		if os.path.exists(path):
			os.remove(path)

	if os.path.exists(Output):
		os.remove(Output)
		
	print 'Do maple stuff'

	os.system(MAPLE_PATH + ' < ' + MplFile)

	template_file = open(TemplateFile, 'r')
	file_write_to = open(Output, 'w')

	print 'Ersetzungen durchfuehren'

	
	try:
		text = template_file.read()
	finally:
		template_file.close()

	i = 0

	for filename in TmpFiles:

		path = pathTmp % filename[0]
		HFile = open(path, 'r')

		text_rep = u'$%d$' % i 

		
		try:
			text_input = HFile.read()
			if filename[1] == 1:
				for str_pair in strings_to_replace:
					text_input = text_input.replace(str_pair[0], str_pair[1])
			
			text = text.replace(text_rep, text_input)
			
		finally:
			HFile.close()

		i = i + 1

	file_write_to.write(text)
	file_write_to.close()

	print 'Temporaere Dateien loeschen'

	for filename in TmpFiles:
		path = pathTmp % filename[0]

		if os.path.exists(path):
			os.remove(path)

if __name__ == "__main__":
	print "Hello World"