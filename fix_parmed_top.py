import subprocess,os,sys

print(sys.argv)

if len(sys.argv)>1:

	input_topology=sys.argv[1]

#	input_topology="topology/system.top"
	topology_file=open(input_topology,"r")
	top_lines=topology_file.readlines()

	out_line=[]
	first_mol=0
	molecule_name=""

	for i in range(0,len(top_lines)):
	    
	    if "[ moleculetype ]" in top_lines[i] or "[ system ]" in top_lines[i]:
	        if molecule_name:
	            if not molecule_name=="SOL" and not molecule_name=="NA" and not molecule_name=="CL":
	            
	                out_line.append("\n#ifdef POSRES\n")
	                out_line.append('#include "posres_%s.itp"\n'%molecule_name)
	                out_line.append("#endif\n\n")

	            if molecule_name=="MOL":
	                out_line.append("\n#ifdef POSRES_LIG\n")
	                out_line.append('#include "posres_%s.itp"\n'%molecule_name)
	                out_line.append("#endif\n\n")

	                
	        
	        first_mol+=1
	        
	        try:
	            posres_file.close()
	        except:
	            pass
	        
	        molecule_name=top_lines[i+2].split()[0]
	                
	        posres_filename='posres_'+molecule_name+".itp"
	        posres_file=open(posres_filename,"w")
	        posres_file.write("[ position_restraints ]\n")
	           
	        
	    if top_lines[i]:
	        if not top_lines[i].startswith(";") and not top_lines[i].startswith("[") and "qtot" in top_lines[i]:
	            arr=top_lines[i].split()
	            if(len(arr))>7:
	                if not float(arr[7])==3:
	                    #print(top_lines[i])
	                    posres_str=str(arr[0])+" 1 2500 2500 2500\n"
	                    posres_file.write(posres_str)
	        
	    
	#    if "[ moleculetype ]" in top_lines[i] and first_mol>1:
	        
	        
	    out_line.append(top_lines[i])
	    
	posres_file.close()


	#print(len(top_lines))

	out_top=open('out_top.top','w')
	for line in out_line:
	    out_top.write(line)
	    
	out_top.close()
else:
	print("Provide a valid topology file\n")
