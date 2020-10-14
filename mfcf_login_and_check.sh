#!/bin/bash

sshpass -p "Mathsprofessor6351*" ssh -tt -o StrictHostKeyChecking=no b2philli@rsubmit.math.private.uwaterloo.ca << EOF
	cd /work/b2philli;
	touch logged_in.txt;
	. submit_jobs.sh;
	exit;
EOF


