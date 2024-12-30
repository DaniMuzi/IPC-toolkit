

current_task=0
while (($current_task < $tasks_to_be_done)); do
	running_tasks=`ps -C IPC --no-headers | wc -l`
	while (($running_tasks < $max_tasks && ${current_task} < ${tasks_to_be_done})); do

		set -- ${_params[$current_task]}
		T=${1}
		N=${2}
		LastFold="T${T}_N${N}"

		n=${3}

		f1="${LastFold}/settings_num${n}.txt"

		f="${full_path}/${f1}"
		# echo "$f"

		./IPC $f &

		running_tasks=`ps -C IPC --no-headers | wc -l`
		((current_task++))

	done
done
wait
