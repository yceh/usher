from multiprocessing.connection import Listener
relevent_sites_file_name = "filtering/data/allRelevantNodesInfSites.txt"
relevent_sites_file_sock = "filtering/data/allRelevantNodesInfSites.sock"
import sys
with Listener(relevent_sites_file_sock, family='AF_UNIX',backlog=int(sys.argv[1])) as listener:
	revelent_sites={}

	with open(relevent_sites_file_name, 'r') as F:
		for line in F:
			line = line.split('\t')
			revelent_sites[(int(line[0]),int(line[1]),int(line[2]))]=line[7].split(',')
	while (True):
		with listener.accept() as conn:
			idx=conn.recv()
			try:
				conn.send(revelent_sites[idx])
			except KeyError:
				print ( "relevent sites"+str(idx)+"not found")
				conn.send([])