import simulate
import train

import sys

def nanosim_simulate():
	simulate.main()

def nanosim_train():
	train.main(sys.argv[1:])
