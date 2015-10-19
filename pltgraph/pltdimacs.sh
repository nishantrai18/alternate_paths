#!/bin/bash

#Works for small graphs in DIMACS format

awk  'BEGIN{
		line=0;
		FS=" ";
		printf  "digraph gitgr {\n";
	}
	{
		if(line<100)
		{
			if($1=="a")
			{
				printf "\t";
				printf "c"$2" ";
				printf "->";
				printf " c"$3;
				printf ";\n";	
			}
		}
		line++;
	}
	END  {
		printf "}\n";
	}' < "rome99.gr" | dot -Tps > out.ps
