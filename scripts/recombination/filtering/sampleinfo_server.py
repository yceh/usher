from multiprocessing.connection import Listener
sampleinfo_file_name = "filtering/data/sampleInfo.txt"
sampleinfo_sock_name = "filtering/data/sampleInfo.sock"
import sys
with Listener(sampleinfo_sock_name, family='AF_UNIX',backlog=int(sys.argv[1])) as listener:
	sample_info={}

	with open(sampleinfo_file_name, 'r') as F:
		for line in F:
			fields=F.split(',')
			samples=fields[1].split(',')
			mutations=fields[2].split(',')
			sample_info[fields[0]]=(samples[:10],mutations)
	while (True):
		with listener.accept() as conn:
			while (True):
				idx=conn.recv()
				if(idx=="END"):
					break
				conn.send(sample_info[idx])