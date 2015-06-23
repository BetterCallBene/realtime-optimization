#!/usr/bin/python

import sys
import GenerateScript as GS

def main(template):
	files = [
			'tmpRTOptFunction', #0
			'tmpRTOptJacobi', #1
			'tmpRTOptHesse' #2
			] 


	if template is None:
		template = 'template_BasisGenQDyn.m'

	file_to = 'BasisGenQDyn.m'

	strings_to_replace = [
						['q(1)', 'q(1, :)'],
						['q(2)', 'q(2, :)'],
						['q(3)', 'q(3, :)'],
						['q(4)', 'q(4, :)'],
						['v(1)', 'v(1, :)'],
						['v(2)', 'v(2, :)'],
						['v(3)', 'v(3, :)'],
						['omega(1)', 'omega(1, :)'],
						['omega(2)', 'omega(2, :)'],
						['omega(3)', 'omega(3, :)'],
						['u(1)', 'u(1, :)'],
						['u(2)', 'u(2, :)'],
						['u(3)', 'u(3, :)'],
						['u(4)', 'u(4, :)'],
						['*', '.*'],
						['/', './'],
						['^', u'.^'],
						['ones', 'onesV']
						]

	file_mpl = 'GenerateJacobiHesse.mpl'

	GS.generate(file_to, files, file_mpl, template, strings_to_replace)

if __name__ == "__main__":
	template = ''
	if sys.argv[1] is not None:
		template = sys.argv[1]

	main(template)