from multiprocessing.connection import Listener
sampleinfo_file_name = "filtering/data/sampleInfo.txt"
sampleinfo_sock_name = "filtering/data/sampleInfo.sock"
import sys
with Listener(sampleinfo_sock_name, family='AF_UNIX',backlog=int(sys.argv[1])) as listener:
	sample_info={}

	with open(sampleinfo_file_name, 'r') as F:
		F.readline()
		for line in F:
			fields=line.split('\t')
			samples=fields[1].split(',')
			mutations=fields[2].split(',')
			sample_info[int(fields[0])]=(samples[:10],mutations)
	while (True):
		with listener.accept() as conn:
			while (True):
				idx=conn.recv()
				if(idx=="END"):
					break
				try:
					conn.send(sample_info[idx])
				except KeyError:
					print("sample_info_server %s requested not found" % idx)
					conn.send(([],[]))