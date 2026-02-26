#!/bin/bash
#SBATCH --job-name LBbe-128
#SBATCH --nodes=128
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:4
#SBATCH --gpus-per-task=1
#SBATCH --time=02:00:00         # Walltime, format: HH:MM:SS
#SBATCH --hint=nomultithread
#SBATCH --mem=0
#SBATCH -A IscrB_MIPLAST
#SBATCH --qos=boost_qos_bprod
#SBATCH -p boost_usr_prod
#SBATCH --exclusive
###SBATCH --qos=boost_qos_dbg
#SBATCH --output=%x_%j.log
#SBATCH --error=stdout
 

#module purge
module load nvhpc/24.3
module load openmpi/4.1.6--nvhpc--24.3

export OMP_NUM_THREADS=1
export OMPI_MCA_io=ompio
export OMPI_MCA_sharedfp=individual



#which mpirun
#mpirun --mca pml ucx --mca btl ^openib --version
#command -v ompi_info &> /dev/null && ompi_info | grep -i ucx || echo "ompi_info non disponibile"

#echo "CUDA_VISIBLE_DEVICES: $CUDA_VISIBLE_DEVICES"

#mycommand="mpirun --mca pml ucx --mca btl ^openib --mca io romio321"

#2>&1

mycommand="mpirun --wdir "$PWD" --bind-to core --map-by ppr:4:node:pe=8 --mca io ompio --mca sharedfp individual "

#$mycommand -np 4 ./main.x 1 1 4 input_single.inp > mio.txt

mynodes="128"  
#mycommand="mpirun --map-by ppr:4:node --mca btl ^openib --mca pml ucx -x UCX_NET_DEVICES='mlx5_0:1' -x CUDA_VISIBLE_DEVICES=0,1,2,3"


if [[ "$mynodes" == "1" ]]; then
    echo "Mode: WEAK - Running weak configuration"
    $mycommand -np 1 ./main_1c_15.x 1 1 1 bench1w.inp > mio1_15_w_1_1_1_1.dat
    $mycommand -np 2 ./main_1c_15.x 1 1 2 bench1w.inp > mio1_15_w_2_1_1_2.dat
    $mycommand -np 4 ./main_1c_15.x 1 2 2 bench1w.inp > mio1_15_w_4_1_2_2.dat
    $mycommand -np 4 ./main_1c_15.x 1 1 4 bench1w.inp > mio1_15_w_4_1_1_4.dat

    $mycommand -np 1 ./main_2c_15.x 1 1 1 benchw.inp > mio2_15_w_1_1_1_1.dat
    $mycommand -np 2 ./main_2c_15.x 1 1 2 benchw.inp > mio2_15_w_2_1_1_2.dat
    $mycommand -np 4 ./main_2c_15.x 1 2 2 benchw.inp > mio2_15_w_4_1_2_2.dat
    $mycommand -np 4 ./main_2c_15.x 1 1 4 benchw.inp > mio2_15_w_4_1_1_4.dat
    
    echo "Mode: STRONG - Running strong configuration"
    
    $mycommand  -np 1 ./main_1c_15.x 1 1 1 bench1.inp > mio1_15_s_1_1_1_1.dat
    $mycommand  -np 2 ./main_1c_15.x 1 1 2 bench1.inp > mio1_15_s_2_1_1_2.dat
    $mycommand  -np 4 ./main_1c_15.x 1 1 4 bench1.inp > mio1_15_s_4_1_1_4.dat
    $mycommand  -np 4 ./main_1c_15.x 1 2 2 bench1.inp > mio1_15_s_4_1_2_2.dat

    $mycommand  -np 1 ./main_2c_15.x 1 1 1 bench.inp > mio2_15_s_1_1_1_1.dat
    $mycommand  -np 2 ./main_2c_15.x 1 1 2 bench.inp > mio2_15_s_2_1_1_2.dat
    $mycommand  -np 4 ./main_2c_15.x 1 1 4 bench.inp > mio2_15_s_4_1_1_4.dat
    $mycommand  -np 4 ./main_2c_15.x 1 2 2 bench.inp > mio2_15_s_4_1_2_2.dat
fi


if [[ "$mynodes" == "2" ]]; then
	echo "Mode: WEAK - Running weak configuration"
        $mycommand -np 8 ./main_1c_15.x 1 1 8 bench1w.inp > mio1_15_w_8_1_1_8.dat
        $mycommand -np 8 ./main_1c_15.x 1 2 4 bench1w.inp > mio1_15_w_8_1_2_4.dat
        $mycommand -np 8 ./main_1c_15.x 2 2 2 bench1w.inp > mio1_15_w_8_2_2_2.dat

	    $mycommand -np 8 ./main_2c_15.x 1 1 8 benchw.inp > mio2_15_w_8_1_1_8.dat
	    $mycommand -np 8 ./main_2c_15.x 1 2 4 benchw.inp > mio2_15_w_8_1_2_4.dat
        $mycommand -np 8 ./main_2c_15.x 2 2 2 benchw.inp > mio2_15_w_8_2_2_2.dat
	    
        echo "Mode: STRONG - Running strong configuration"
        $mycommand -np 8 ./main_1c_15.x  2 2 2 bench1.inp > mio1_15_s_8_2_2_2.dat
        $mycommand -np 8 ./main_1c_15.x  1 2 4 bench1.inp > mio1_15_s_8_1_2_4.dat

        $mycommand -np 8 ./main_2c_15.x  2 2 2 bench.inp > mio2_15_s_8_2_2_2.dat
        $mycommand -np 8 ./main_2c_15.x  1 2 4 bench.inp > mio2_15_s_8_1_2_4.dat

fi

if [[ "$mynodes" == "4" ]]; then
        echo "Mode: WEAK - Running weak configuration"

        $mycommand -np 16 ./main_1c_15.x 1 1 16 bench1w.inp > mio1_15_w_16_1_1_16.dat
        $mycommand -np 16 ./main_1c_15.x 1 4 4 bench1w.inp > mio1_15_w_16_1_4_4.dat
        $mycommand -np 16 ./main_1c_15.x 2 2 4 bench1w.inp > mio1_15_w_16_2_2_4.dat

	    $mycommand -np 16 ./main_2c_15.x 1 1 16 benchw.inp > mio2_15_w_16_1_1_16.dat
	    $mycommand -np 16 ./main_2c_15.x 1 4 4 benchw.inp > mio2_15_w_16_1_4_4.dat
        $mycommand -np 16 ./main_2c_15.x 2 2 4 benchw.inp > mio2_15_w_16_2_2_4.dat
	    
        echo "Mode: STRONG - Running strong configuration"
       $mycommand -np 16 ./main_1c_15.x 2 2 4 bench1.inp > mio1_15_s_16_2_2_4.dat
       $mycommand -np 16 ./main_1c_15.x 1 4 4 bench1.inp > mio1_15_s_16_1_4_4.dat

       $mycommand -np 16 ./main_2c_15.x 2 2 4 bench.inp > mio2_15_s_16_2_2_4.dat
       $mycommand -np 16 ./main_2c_15.x 1 4 4 bench.inp > mio2_15_s_16_1_4_4.dat

fi

if [[ "$mynodes" == "8" ]]; then

    echo "Mode: WEAK - Running weak configuration"
    
    $mycommand -np 32 ./main_1c_15.x 1 1 32 bench1w.inp > mio1_15_w_32_1_1_32.dat
    $mycommand -np 32 ./main_1c_15.x 1 4 8 bench1w.inp > mio1_15_w_32_1_4_8.dat
    $mycommand -np 32 ./main_1c_15.x 2 4 4 bench1w.inp > mio1_15_w_32_2_4_4.dat
    
    $mycommand -np 32 ./main_2c_15.x 1 1 32 benchw.inp > mio2_15_w_32_1_1_32.dat
    $mycommand -np 32 ./main_2c_15.x 1 4 8 benchw.inp > mio2_15_w_32_1_4_8.dat
    $mycommand -np 32 ./main_2c_15.x 2 4 4 benchw.inp > mio2_15_w_32_2_4_4.dat
    
    echo "Mode: STRONG - Running strong configuration"
    $mycommand -np 32 ./main_1c_15.x 2 4 4 bench1.inp > mio1_15_s_32_2_4_4.dat
    $mycommand -np 32 ./main_1c_15.x 1 4 8 bench1.inp > mio1_15_s_32_1_4_8.dat
    
    $mycommand -np 32 ./main_2c_15.x 2 4 4 bench.inp > mio2_15_s_32_2_4_4.dat
    $mycommand -np 32 ./main_2c_15.x 1 4 8 bench.inp > mio2_15_s_32_1_4_8.dat

fi

if [[ "$mynodes" == "16" ]]; then
    echo "Mode: WEAK - Running weak configuration"
    $mycommand -np 64 ./main_1c_15.x  1 1 64 bench1w.inp > mio1_15_w_64_1_1_64.dat
    $mycommand -np 64 ./main_1c_15.x  1 8 8 bench1w.inp > mio1_15_w_64_1_8_8.dat
    $mycommand -np 64 ./main_1c_15.x  2 4 8 bench1w.inp > mio1_15_w_64_2_4_8.dat
    $mycommand -np 64 ./main_1c_15.x  4 4 4 bench1w.inp > mio1_15_w_64_4_4_4.dat
    
    $mycommand -np 64 ./main_2c_15.x  1 1 64 benchw.inp > mio2_15_w_64_1_1_64.dat
    $mycommand -np 64 ./main_2c_15.x  1 8 8 benchw.inp > mio2_15_w_64_1_8_8.dat
    $mycommand -np 64 ./main_2c_15.x  2 4 8 benchw.inp > mio2_15_w_64_2_4_8.dat
    $mycommand -np 64 ./main_2c_15.x  4 4 4 benchw.inp > mio2_15_w_64_4_4_4.dat
    
    echo "Mode: STRONG - Running strong configuration"
    $mycommand -np 64 ./main_1c_15.x 4 4 4 bench1.inp > mio1_15_s_64_4_4_4.dat
    $mycommand -np 64 ./main_1c_15.x 2 4 8 bench1.inp > mio1_15_s_64_2_4_8.dat
    $mycommand -np 64 ./main_1c_15.x 1 8 8 bench1.inp > mio1_15_s_64_1_8_8.dat
    
    $mycommand -np 64 ./main_2c_15.x 4 4 4 bench.inp > mio2_15_s_64_4_4_4.dat
    $mycommand -np 64 ./main_2c_15.x 2 4 8 bench.inp > mio2_15_s_64_2_4_8.dat
    $mycommand -np 64 ./main_2c_15.x 1 8 8 bench.inp > mio2_15_s_64_1_8_8.dat

fi

if [[ "$mynodes" == "1" ]]; then
    echo "Mode: WEAK - Running weak configuration"
    $mycommand -np 1 ./main_1c_19.x 1 1 1 bench1w.inp > mio1_19_w_1_1_1_1.dat
    $mycommand -np 2 ./main_1c_19.x 1 1 2 bench1w.inp > mio1_19_w_2_1_1_2.dat
    $mycommand -np 4 ./main_1c_19.x 1 2 2 bench1w.inp > mio1_19_w_4_1_2_2.dat
    $mycommand -np 4 ./main_1c_19.x 1 1 4 bench1w.inp > mio1_19_w_4_1_1_4.dat

    $mycommand -np 1 ./main_2c_19.x 1 1 1 benchw.inp > mio2_19_w_1_1_1_1.dat
    $mycommand -np 2 ./main_2c_19.x 1 1 2 benchw.inp > mio2_19_w_2_1_1_2.dat
    $mycommand -np 4 ./main_2c_19.x 1 2 2 benchw.inp > mio2_19_w_4_1_2_2.dat
    $mycommand -np 4 ./main_2c_19.x 1 1 4 benchw.inp > mio2_19_w_4_1_1_4.dat
    
    echo "Mode: STRONG - Running strong configuration"
    
    $mycommand  -np 1 ./main_1c_19.x 1 1 1 bench1.inp > mio1_19_s_1_1_1_1.dat
    $mycommand  -np 2 ./main_1c_19.x 1 1 2 bench1.inp > mio1_19_s_2_1_1_2.dat
    $mycommand  -np 4 ./main_1c_19.x 1 1 4 bench1.inp > mio1_19_s_4_1_1_4.dat
    $mycommand  -np 4 ./main_1c_19.x 1 2 2 bench1.inp > mio1_19_s_4_1_2_2.dat

    $mycommand  -np 1 ./main_2c_19.x 1 1 1 bench.inp > mio2_19_s_1_1_1_1.dat
    $mycommand  -np 2 ./main_2c_19.x 1 1 2 bench.inp > mio2_19_s_2_1_1_2.dat
    $mycommand  -np 4 ./main_2c_19.x 1 1 4 bench.inp > mio2_19_s_4_1_1_4.dat
    $mycommand  -np 4 ./main_2c_19.x 1 2 2 bench.inp > mio2_19_s_4_1_2_2.dat
fi


if [[ "$mynodes" == "2" ]]; then
	echo "Mode: WEAK - Running weak configuration"
        $mycommand -np 8 ./main_1c_19.x 1 1 8 bench1w.inp > mio1_19_w_8_1_1_8.dat
        $mycommand -np 8 ./main_1c_19.x 1 2 4 bench1w.inp > mio1_19_w_8_1_2_4.dat
        $mycommand -np 8 ./main_1c_19.x 2 2 2 bench1w.inp > mio1_19_w_8_2_2_2.dat

	    $mycommand -np 8 ./main_2c_19.x 1 1 8 benchw.inp > mio2_19_w_8_1_1_8.dat
	    $mycommand -np 8 ./main_2c_19.x 1 2 4 benchw.inp > mio2_19_w_8_1_2_4.dat
        $mycommand -np 8 ./main_2c_19.x 2 2 2 benchw.inp > mio2_19_w_8_2_2_2.dat
	    
        echo "Mode: STRONG - Running strong configuration"
        $mycommand -np 8 ./main_1c_19.x  2 2 2 bench1.inp > mio1_19_s_8_2_2_2.dat
        $mycommand -np 8 ./main_1c_19.x  1 2 4 bench1.inp > mio1_19_s_8_1_2_4.dat

        $mycommand -np 8 ./main_2c_19.x  2 2 2 bench.inp > mio2_19_s_8_2_2_2.dat
        $mycommand -np 8 ./main_2c_19.x  1 2 4 bench.inp > mio2_19_s_8_1_2_4.dat

fi

if [[ "$mynodes" == "4" ]]; then
        echo "Mode: WEAK - Running weak configuration"

        $mycommand -np 16 ./main_1c_19.x 1 1 16 bench1w.inp > mio1_19_w_16_1_1_16.dat
        $mycommand -np 16 ./main_1c_19.x 1 4 4 bench1w.inp > mio1_19_w_16_1_4_4.dat
        $mycommand -np 16 ./main_1c_19.x 2 2 4 bench1w.inp > mio1_19_w_16_2_2_4.dat

	    $mycommand -np 16 ./main_2c_19.x 1 1 16 benchw.inp > mio2_19_w_16_1_1_16.dat
	    $mycommand -np 16 ./main_2c_19.x 1 4 4 benchw.inp > mio2_19_w_16_1_4_4.dat
        $mycommand -np 16 ./main_2c_19.x 2 2 4 benchw.inp > mio2_19_w_16_2_2_4.dat
	    
        echo "Mode: STRONG - Running strong configuration"
       $mycommand -np 16 ./main_1c_19.x 2 2 4 bench1.inp > mio1_19_s_16_2_2_4.dat
       $mycommand -np 16 ./main_1c_19.x 1 4 4 bench1.inp > mio1_19_s_16_1_4_4.dat

       $mycommand -np 16 ./main_2c_19.x 2 2 4 bench.inp > mio2_19_s_16_2_2_4.dat
       $mycommand -np 16 ./main_2c_19.x 1 4 4 bench.inp > mio2_19_s_16_1_4_4.dat

fi

if [[ "$mynodes" == "8" ]]; then

    echo "Mode: WEAK - Running weak configuration"
    
    $mycommand -np 32 ./main_1c_19.x 1 1 32 bench1w.inp > mio1_19_w_32_1_1_32.dat
    $mycommand -np 32 ./main_1c_19.x 1 4 8 bench1w.inp > mio1_19_w_32_1_4_8.dat
    $mycommand -np 32 ./main_1c_19.x 2 4 4 bench1w.inp > mio1_19_w_32_2_4_4.dat
    
    $mycommand -np 32 ./main_2c_19.x 1 1 32 benchw.inp > mio2_19_w_32_1_1_32.dat
    $mycommand -np 32 ./main_2c_19.x 1 4 8 benchw.inp > mio2_19_w_32_1_4_8.dat
    $mycommand -np 32 ./main_2c_19.x 2 4 4 benchw.inp > mio2_19_w_32_2_4_4.dat
    
    echo "Mode: STRONG - Running strong configuration"
    $mycommand -np 32 ./main_1c_19.x 2 4 4 bench1.inp > mio1_19_s_32_2_4_4.dat
    $mycommand -np 32 ./main_1c_19.x 1 4 8 bench1.inp > mio1_19_s_32_1_4_8.dat
    
    $mycommand -np 32 ./main_2c_19.x 2 4 4 bench.inp > mio2_19_s_32_2_4_4.dat
    $mycommand -np 32 ./main_2c_19.x 1 4 8 bench.inp > mio2_19_s_32_1_4_8.dat

fi

if [[ "$mynodes" == "16" ]]; then
    echo "Mode: WEAK - Running weak configuration"
    $mycommand -np 64 ./main_1c_19.x  1 1 64 bench1w.inp > mio1_19_w_64_1_1_64.dat
    $mycommand -np 64 ./main_1c_19.x  1 8 8 bench1w.inp > mio1_19_w_64_1_8_8.dat
    $mycommand -np 64 ./main_1c_19.x  2 4 8 bench1w.inp > mio1_19_w_64_2_4_8.dat
    $mycommand -np 64 ./main_1c_19.x  4 4 4 bench1w.inp > mio1_19_w_64_4_4_4.dat
    
    $mycommand -np 64 ./main_2c_19.x  1 1 64 benchw.inp > mio2_19_w_64_1_1_64.dat
    $mycommand -np 64 ./main_2c_19.x  1 8 8 benchw.inp > mio2_19_w_64_1_8_8.dat
    $mycommand -np 64 ./main_2c_19.x  2 4 8 benchw.inp > mio2_19_w_64_2_4_8.dat
    $mycommand -np 64 ./main_2c_19.x  4 4 4 benchw.inp > mio2_19_w_64_4_4_4.dat
    
    echo "Mode: STRONG - Running strong configuration"
    $mycommand -np 64 ./main_1c_19.x 4 4 4 bench1.inp > mio1_19_s_64_4_4_4.dat
    $mycommand -np 64 ./main_1c_19.x 2 4 8 bench1.inp > mio1_19_s_64_2_4_8.dat
    $mycommand -np 64 ./main_1c_19.x 1 8 8 bench1.inp > mio1_19_s_64_1_8_8.dat
    
    $mycommand -np 64 ./main_2c_19.x 4 4 4 bench.inp > mio2_19_s_64_4_4_4.dat
    $mycommand -np 64 ./main_2c_19.x 2 4 8 bench.inp > mio2_19_s_64_2_4_8.dat
    $mycommand -np 64 ./main_2c_19.x 1 8 8 bench.inp > mio2_19_s_64_1_8_8.dat

fi


if [[ "$mynodes" == "1" ]]; then
    echo "Mode: WEAK - Running weak configuration"
    $mycommand -np 1 ./main_1c_27.x 1 1 1 bench1w.inp > mio1_27_w_1_1_1_1.dat
    $mycommand -np 2 ./main_1c_27.x 1 1 2 bench1w.inp > mio1_27_w_2_1_1_2.dat
    $mycommand -np 4 ./main_1c_27.x 1 2 2 bench1w.inp > mio1_27_w_4_1_2_2.dat
    $mycommand -np 4 ./main_1c_27.x 1 1 4 bench1w.inp > mio1_27_w_4_1_1_4.dat

    $mycommand -np 1 ./main_2c_27.x 1 1 1 benchw.inp > mio2_27_w_1_1_1_1.dat
    $mycommand -np 2 ./main_2c_27.x 1 1 2 benchw.inp > mio2_27_w_2_1_1_2.dat
    $mycommand -np 4 ./main_2c_27.x 1 2 2 benchw.inp > mio2_27_w_4_1_2_2.dat
    $mycommand -np 4 ./main_2c_27.x 1 1 4 benchw.inp > mio2_27_w_4_1_1_4.dat
    
    echo "Mode: STRONG - Running strong configuration"
    
    $mycommand  -np 1 ./main_1c_27.x 1 1 1 bench1.inp > mio1_27_s_1_1_1_1.dat
    $mycommand  -np 2 ./main_1c_27.x 1 1 2 bench1.inp > mio1_27_s_2_1_1_2.dat
    $mycommand  -np 4 ./main_1c_27.x 1 1 4 bench1.inp > mio1_27_s_4_1_1_4.dat
    $mycommand  -np 4 ./main_1c_27.x 1 2 2 bench1.inp > mio1_27_s_4_1_2_2.dat

    $mycommand  -np 1 ./main_2c_27.x 1 1 1 bench.inp > mio2_27_s_1_1_1_1.dat
    $mycommand  -np 2 ./main_2c_27.x 1 1 2 bench.inp > mio2_27_s_2_1_1_2.dat
    $mycommand  -np 4 ./main_2c_27.x 1 1 4 bench.inp > mio2_27_s_4_1_1_4.dat
    $mycommand  -np 4 ./main_2c_27.x 1 2 2 bench.inp > mio2_27_s_4_1_2_2.dat
fi


if [[ "$mynodes" == "2" ]]; then
	echo "Mode: WEAK - Running weak configuration"
        $mycommand -np 8 ./main_1c_27.x 1 1 8 bench1w.inp > mio1_27_w_8_1_1_8.dat
        $mycommand -np 8 ./main_1c_27.x 1 2 4 bench1w.inp > mio1_27_w_8_1_2_4.dat
        $mycommand -np 8 ./main_1c_27.x 2 2 2 bench1w.inp > mio1_27_w_8_2_2_2.dat

	    $mycommand -np 8 ./main_2c_27.x 1 1 8 benchw.inp > mio2_27_w_8_1_1_8.dat
	    $mycommand -np 8 ./main_2c_27.x 1 2 4 benchw.inp > mio2_27_w_8_1_2_4.dat
        $mycommand -np 8 ./main_2c_27.x 2 2 2 benchw.inp > mio2_27_w_8_2_2_2.dat
	    
        echo "Mode: STRONG - Running strong configuration"
        $mycommand -np 8 ./main_1c_27.x  2 2 2 bench1.inp > mio1_27_s_8_2_2_2.dat
        $mycommand -np 8 ./main_1c_27.x  1 2 4 bench1.inp > mio1_27_s_8_1_2_4.dat

        $mycommand -np 8 ./main_2c_27.x  2 2 2 bench.inp > mio2_27_s_8_2_2_2.dat
        $mycommand -np 8 ./main_2c_27.x  1 2 4 bench.inp > mio2_27_s_8_1_2_4.dat

fi

if [[ "$mynodes" == "4" ]]; then
        echo "Mode: WEAK - Running weak configuration"

        $mycommand -np 16 ./main_1c_27.x 1 1 16 bench1w.inp > mio1_27_w_16_1_1_16.dat
        $mycommand -np 16 ./main_1c_27.x 1 4 4 bench1w.inp > mio1_27_w_16_1_4_4.dat
        $mycommand -np 16 ./main_1c_27.x 2 2 4 bench1w.inp > mio1_27_w_16_2_2_4.dat

	    $mycommand -np 16 ./main_2c_27.x 1 1 16 benchw.inp > mio2_27_w_16_1_1_16.dat
	    $mycommand -np 16 ./main_2c_27.x 1 4 4 benchw.inp > mio2_27_w_16_1_4_4.dat
        $mycommand -np 16 ./main_2c_27.x 2 2 4 benchw.inp > mio2_27_w_16_2_2_4.dat
	    
        echo "Mode: STRONG - Running strong configuration"
       $mycommand -np 16 ./main_1c_27.x 2 2 4 bench1.inp > mio1_27_s_16_2_2_4.dat
       $mycommand -np 16 ./main_1c_27.x 1 4 4 bench1.inp > mio1_27_s_16_1_4_4.dat

       $mycommand -np 16 ./main_2c_27.x 2 2 4 bench.inp > mio2_27_s_16_2_2_4.dat
       $mycommand -np 16 ./main_2c_27.x 1 4 4 bench.inp > mio2_27_s_16_1_4_4.dat

fi

if [[ "$mynodes" == "8" ]]; then

    echo "Mode: WEAK - Running weak configuration"
    
    $mycommand -np 32 ./main_1c_27.x 1 1 32 bench1w.inp > mio1_27_w_32_1_1_32.dat
    $mycommand -np 32 ./main_1c_27.x 1 4 8 bench1w.inp > mio1_27_w_32_1_4_8.dat
    $mycommand -np 32 ./main_1c_27.x 2 4 4 bench1w.inp > mio1_27_w_32_2_4_4.dat
    
    $mycommand -np 32 ./main_2c_27.x 1 1 32 benchw.inp > mio2_27_w_32_1_1_32.dat
    $mycommand -np 32 ./main_2c_27.x 1 4 8 benchw.inp > mio2_27_w_32_1_4_8.dat
    $mycommand -np 32 ./main_2c_27.x 2 4 4 benchw.inp > mio2_27_w_32_2_4_4.dat
    
    echo "Mode: STRONG - Running strong configuration"
    $mycommand -np 32 ./main_1c_27.x 2 4 4 bench1.inp > mio1_27_s_32_2_4_4.dat
    $mycommand -np 32 ./main_1c_27.x 1 4 8 bench1.inp > mio1_27_s_32_1_4_8.dat
    
    $mycommand -np 32 ./main_2c_27.x 2 4 4 bench.inp > mio2_27_s_32_2_4_4.dat
    $mycommand -np 32 ./main_2c_27.x 1 4 8 bench.inp > mio2_27_s_32_1_4_8.dat

fi

if [[ "$mynodes" == "16" ]]; then
    echo "Mode: WEAK - Running weak configuration"
    $mycommand -np 64 ./main_1c_27.x  1 1 64 bench1w.inp > mio1_27_w_64_1_1_64.dat
    $mycommand -np 64 ./main_1c_27.x  1 8 8 bench1w.inp > mio1_27_w_64_1_8_8.dat
    $mycommand -np 64 ./main_1c_27.x  2 4 8 bench1w.inp > mio1_27_w_64_2_4_8.dat
    $mycommand -np 64 ./main_1c_27.x  4 4 4 bench1w.inp > mio1_27_w_64_4_4_4.dat
    
    $mycommand -np 64 ./main_2c_27.x  1 1 64 benchw.inp > mio2_27_w_64_1_1_64.dat
    $mycommand -np 64 ./main_2c_27.x  1 8 8 benchw.inp > mio2_27_w_64_1_8_8.dat
    $mycommand -np 64 ./main_2c_27.x  2 4 8 benchw.inp > mio2_27_w_64_2_4_8.dat
    $mycommand -np 64 ./main_2c_27.x  4 4 4 benchw.inp > mio2_27_w_64_4_4_4.dat
    
    echo "Mode: STRONG - Running strong configuration"
    $mycommand -np 64 ./main_1c_27.x 4 4 4 bench1.inp > mio1_27_s_64_4_4_4.dat
    $mycommand -np 64 ./main_1c_27.x 2 4 8 bench1.inp > mio1_27_s_64_2_4_8.dat
    $mycommand -np 64 ./main_1c_27.x 1 8 8 bench1.inp > mio1_27_s_64_1_8_8.dat
    
    $mycommand -np 64 ./main_2c_27.x 4 4 4 bench.inp > mio2_27_s_64_4_4_4.dat
    $mycommand -np 64 ./main_2c_27.x 2 4 8 bench.inp > mio2_27_s_64_2_4_8.dat
    $mycommand -np 64 ./main_2c_27.x 1 8 8 bench.inp > mio2_27_s_64_1_8_8.dat

fi


if [[ "$mynodes" == "1" ]]; then
    echo "Mode: WEAK - Running weak configuration"
    $mycommand -np 1 ./main_1c_27high.x 1 1 1 bench1w.inp > mio1_27h_w_1_1_1_1.dat
    $mycommand -np 2 ./main_1c_27high.x 1 1 2 bench1w.inp > mio1_27h_w_2_1_1_2.dat
    $mycommand -np 4 ./main_1c_27high.x 1 2 2 bench1w.inp > mio1_27h_w_4_1_2_2.dat
    $mycommand -np 4 ./main_1c_27high.x 1 1 4 bench1w.inp > mio1_27h_w_4_1_1_4.dat

    $mycommand -np 1 ./main_2c_27high.x 1 1 1 benchw.inp > mio2_27h_w_1_1_1_1.dat
    $mycommand -np 2 ./main_2c_27high.x 1 1 2 benchw.inp > mio2_27h_w_2_1_1_2.dat
    $mycommand -np 4 ./main_2c_27high.x 1 2 2 benchw.inp > mio2_27h_w_4_1_2_2.dat
    $mycommand -np 4 ./main_2c_27high.x 1 1 4 benchw.inp > mio2_27h_w_4_1_1_4.dat
    
    echo "Mode: STRONG - Running strong configuration"
    
    $mycommand  -np 1 ./main_1c_27high.x 1 1 1 bench1.inp > mio1_27h_s_1_1_1_1.dat
    $mycommand  -np 2 ./main_1c_27high.x 1 1 2 bench1.inp > mio1_27h_s_2_1_1_2.dat
    $mycommand  -np 4 ./main_1c_27high.x 1 1 4 bench1.inp > mio1_27h_s_4_1_1_4.dat
    $mycommand  -np 4 ./main_1c_27high.x 1 2 2 bench1.inp > mio1_27h_s_4_1_2_2.dat

    $mycommand  -np 1 ./main_2c_27high.x 1 1 1 bench.inp > mio2_27h_s_1_1_1_1.dat
    $mycommand  -np 2 ./main_2c_27high.x 1 1 2 bench.inp > mio2_27h_s_2_1_1_2.dat
    $mycommand  -np 4 ./main_2c_27high.x 1 1 4 bench.inp > mio2_27h_s_4_1_1_4.dat
    $mycommand  -np 4 ./main_2c_27high.x 1 2 2 bench.inp > mio2_27h_s_4_1_2_2.dat
fi


if [[ "$mynodes" == "2" ]]; then
	echo "Mode: WEAK - Running weak configuration"
        $mycommand -np 8 ./main_1c_27high.x 1 1 8 bench1w.inp > mio1_27h_w_8_1_1_8.dat
        $mycommand -np 8 ./main_1c_27high.x 1 2 4 bench1w.inp > mio1_27h_w_8_1_2_4.dat
        $mycommand -np 8 ./main_1c_27high.x 2 2 2 bench1w.inp > mio1_27h_w_8_2_2_2.dat

	    $mycommand -np 8 ./main_2c_27high.x 1 1 8 benchw.inp > mio2_27h_w_8_1_1_8.dat
	    $mycommand -np 8 ./main_2c_27high.x 1 2 4 benchw.inp > mio2_27h_w_8_1_2_4.dat
        $mycommand -np 8 ./main_2c_27high.x 2 2 2 benchw.inp > mio2_27h_w_8_2_2_2.dat
	    
        echo "Mode: STRONG - Running strong configuration"
        $mycommand -np 8 ./main_1c_27high.x  2 2 2 bench1.inp > mio1_27h_s_8_2_2_2.dat
        $mycommand -np 8 ./main_1c_27high.x  1 2 4 bench1.inp > mio1_27h_s_8_1_2_4.dat

        $mycommand -np 8 ./main_2c_27high.x  2 2 2 bench.inp > mio2_27h_s_8_2_2_2.dat
        $mycommand -np 8 ./main_2c_27high.x  1 2 4 bench.inp > mio2_27h_s_8_1_2_4.dat

fi

if [[ "$mynodes" == "4" ]]; then
        echo "Mode: WEAK - Running weak configuration"

        $mycommand -np 16 ./main_1c_27high.x 1 1 16 bench1w.inp > mio1_27h_w_16_1_1_16.dat
        $mycommand -np 16 ./main_1c_27high.x 1 4 4 bench1w.inp > mio1_27h_w_16_1_4_4.dat
        $mycommand -np 16 ./main_1c_27high.x 2 2 4 bench1w.inp > mio1_27h_w_16_2_2_4.dat

	    $mycommand -np 16 ./main_2c_27high.x 1 1 16 benchw.inp > mio2_27h_w_16_1_1_16.dat
	    $mycommand -np 16 ./main_2c_27high.x 1 4 4 benchw.inp > mio2_27h_w_16_1_4_4.dat
        $mycommand -np 16 ./main_2c_27high.x 2 2 4 benchw.inp > mio2_27h_w_16_2_2_4.dat
	    
        echo "Mode: STRONG - Running strong configuration"
       $mycommand -np 16 ./main_1c_27high.x 2 2 4 bench1.inp > mio1_27h_s_16_2_2_4.dat
       $mycommand -np 16 ./main_1c_27high.x 1 4 4 bench1.inp > mio1_27h_s_16_1_4_4.dat

       $mycommand -np 16 ./main_2c_27high.x 2 2 4 bench.inp > mio2_27h_s_16_2_2_4.dat
       $mycommand -np 16 ./main_2c_27high.x 1 4 4 bench.inp > mio2_27h_s_16_1_4_4.dat

fi

if [[ "$mynodes" == "8" ]]; then

    echo "Mode: WEAK - Running weak configuration"
    
    $mycommand -np 32 ./main_1c_27high.x 1 1 32 bench1w.inp > mio1_27h_w_32_1_1_32.dat
    $mycommand -np 32 ./main_1c_27high.x 1 4 8 bench1w.inp > mio1_27h_w_32_1_4_8.dat
    $mycommand -np 32 ./main_1c_27high.x 2 4 4 bench1w.inp > mio1_27h_w_32_2_4_4.dat
    
    $mycommand -np 32 ./main_2c_27high.x 1 1 32 benchw.inp > mio2_27h_w_32_1_1_32.dat
    $mycommand -np 32 ./main_2c_27high.x 1 4 8 benchw.inp > mio2_27h_w_32_1_4_8.dat
    $mycommand -np 32 ./main_2c_27high.x 2 4 4 benchw.inp > mio2_27h_w_32_2_4_4.dat
    
    echo "Mode: STRONG - Running strong configuration"
    $mycommand -np 32 ./main_1c_27high.x 2 4 4 bench1.inp > mio1_27h_s_32_2_4_4.dat
    $mycommand -np 32 ./main_1c_27high.x 1 4 8 bench1.inp > mio1_27h_s_32_1_4_8.dat
    
    $mycommand -np 32 ./main_2c_27high.x 2 4 4 bench.inp > mio2_27h_s_32_2_4_4.dat
    $mycommand -np 32 ./main_2c_27high.x 1 4 8 bench.inp > mio2_27h_s_32_1_4_8.dat

fi

if [[ "$mynodes" == "16" ]]; then
    echo "Mode: WEAK - Running weak configuration"
    $mycommand -np 64 ./main_1c_27high.x  1 1 64 bench1w.inp > mio1_27h_w_64_1_1_64.dat
    $mycommand -np 64 ./main_1c_27high.x  1 8 8 bench1w.inp > mio1_27h_w_64_1_8_8.dat
    $mycommand -np 64 ./main_1c_27high.x  2 4 8 bench1w.inp > mio1_27h_w_64_2_4_8.dat
    $mycommand -np 64 ./main_1c_27high.x  4 4 4 bench1w.inp > mio1_27h_w_64_4_4_4.dat
    
    $mycommand -np 64 ./main_2c_27high.x  1 1 64 benchw.inp > mio2_27h_w_64_1_1_64.dat
    $mycommand -np 64 ./main_2c_27high.x  1 8 8 benchw.inp > mio2_27h_w_64_1_8_8.dat
    $mycommand -np 64 ./main_2c_27high.x  2 4 8 benchw.inp > mio2_27h_w_64_2_4_8.dat
    $mycommand -np 64 ./main_2c_27high.x  4 4 4 benchw.inp > mio2_27h_w_64_4_4_4.dat
    
    echo "Mode: STRONG - Running strong configuration"
    $mycommand -np 64 ./main_1c_27high.x 4 4 4 bench1.inp > mio1_27h_s_64_4_4_4.dat
    $mycommand -np 64 ./main_1c_27high.x 2 4 8 bench1.inp > mio1_27h_s_64_2_4_8.dat
    $mycommand -np 64 ./main_1c_27high.x 1 8 8 bench1.inp > mio1_27h_s_64_1_8_8.dat
    
    $mycommand -np 64 ./main_2c_27high.x 4 4 4 bench.inp > mio2_27h_s_64_4_4_4.dat
    $mycommand -np 64 ./main_2c_27high.x 2 4 8 bench.inp > mio2_27h_s_64_2_4_8.dat
    $mycommand -np 64 ./main_2c_27high.x 1 8 8 bench.inp > mio2_27h_s_64_1_8_8.dat

fi

###########################


if [[ "$mynodes" == "1" ]]; then
    echo "Mode: WEAK - Running weak configuration"
    $mycommand -np 1 ./main_1c_15_d.x 1 1 1 bench1w.inp > mio1_15_w_1_1_1_1_d.dat
    $mycommand -np 2 ./main_1c_15_d.x 1 1 2 bench1w.inp > mio1_15_w_2_1_1_2_d.dat
    $mycommand -np 4 ./main_1c_15_d.x 1 2 2 bench1w.inp > mio1_15_w_4_1_2_2_d.dat
    $mycommand -np 4 ./main_1c_15_d.x 1 1 4 bench1w.inp > mio1_15_w_4_1_1_4_d.dat

    $mycommand -np 1 ./main_2c_15_d.x 1 1 1 benchw.inp > mio2_15_w_1_1_1_1_d.dat
    $mycommand -np 2 ./main_2c_15_d.x 1 1 2 benchw.inp > mio2_15_w_2_1_1_2_d.dat
    $mycommand -np 4 ./main_2c_15_d.x 1 2 2 benchw.inp > mio2_15_w_4_1_2_2_d.dat
    $mycommand -np 4 ./main_2c_15_d.x 1 1 4 benchw.inp > mio2_15_w_4_1_1_4_d.dat
    
    echo "Mode: STRONG - Running strong configuration"
    
    $mycommand  -np 1 ./main_1c_15_d.x 1 1 1 bench1.inp > mio1_15_s_1_1_1_1_d.dat
    $mycommand  -np 2 ./main_1c_15_d.x 1 1 2 bench1.inp > mio1_15_s_2_1_1_2_d.dat
    $mycommand  -np 4 ./main_1c_15_d.x 1 1 4 bench1.inp > mio1_15_s_4_1_1_4_d.dat
    $mycommand  -np 4 ./main_1c_15_d.x 1 2 2 bench1.inp > mio1_15_s_4_1_2_2_d.dat

    $mycommand  -np 1 ./main_2c_15_d.x 1 1 1 bench.inp > mio2_15_s_1_1_1_1_d.dat
    $mycommand  -np 2 ./main_2c_15_d.x 1 1 2 bench.inp > mio2_15_s_2_1_1_2_d.dat
    $mycommand  -np 4 ./main_2c_15_d.x 1 1 4 bench.inp > mio2_15_s_4_1_1_4_d.dat
    $mycommand  -np 4 ./main_2c_15_d.x 1 2 2 bench.inp > mio2_15_s_4_1_2_2_d.dat
fi


if [[ "$mynodes" == "2" ]]; then
	echo "Mode: WEAK - Running weak configuration"
        $mycommand -np 8 ./main_1c_15_d.x 1 1 8 bench1w.inp > mio1_15_w_8_1_1_8_d.dat
        $mycommand -np 8 ./main_1c_15_d.x 1 2 4 bench1w.inp > mio1_15_w_8_1_2_4_d.dat
        $mycommand -np 8 ./main_1c_15_d.x 2 2 2 bench1w.inp > mio1_15_w_8_2_2_2_d.dat

	    $mycommand -np 8 ./main_2c_15_d.x 1 1 8 benchw.inp > mio2_15_w_8_1_1_8_d.dat
	    $mycommand -np 8 ./main_2c_15_d.x 1 2 4 benchw.inp > mio2_15_w_8_1_2_4_d.dat
        $mycommand -np 8 ./main_2c_15_d.x 2 2 2 benchw.inp > mio2_15_w_8_2_2_2_d.dat
	    
        echo "Mode: STRONG - Running strong configuration"
        $mycommand -np 8 ./main_1c_15_d.x  2 2 2 bench1.inp > mio1_15_s_8_2_2_2_d.dat
        $mycommand -np 8 ./main_1c_15_d.x  1 2 4 bench1.inp > mio1_15_s_8_1_2_4_d.dat

        $mycommand -np 8 ./main_2c_15_d.x  2 2 2 bench.inp > mio2_15_s_8_2_2_2_d.dat
        $mycommand -np 8 ./main_2c_15_d.x  1 2 4 bench.inp > mio2_15_s_8_1_2_4_d.dat

fi

if [[ "$mynodes" == "4" ]]; then
        echo "Mode: WEAK - Running weak configuration"

        $mycommand -np 16 ./main_1c_15_d.x 1 1 16 bench1w.inp > mio1_15_w_16_1_1_16_d.dat
        $mycommand -np 16 ./main_1c_15_d.x 1 4 4 bench1w.inp > mio1_15_w_16_1_4_4_d.dat
        $mycommand -np 16 ./main_1c_15_d.x 2 2 4 bench1w.inp > mio1_15_w_16_2_2_4_d.dat

	    $mycommand -np 16 ./main_2c_15_d.x 1 1 16 benchw.inp > mio2_15_w_16_1_1_16_d.dat
	    $mycommand -np 16 ./main_2c_15_d.x 1 4 4 benchw.inp > mio2_15_w_16_1_4_4_d.dat
        $mycommand -np 16 ./main_2c_15_d.x 2 2 4 benchw.inp > mio2_15_w_16_2_2_4_d.dat
	    
        echo "Mode: STRONG - Running strong configuration"
       $mycommand -np 16 ./main_1c_15_d.x 2 2 4 bench1.inp > mio1_15_s_16_2_2_4_d.dat
       $mycommand -np 16 ./main_1c_15_d.x 1 4 4 bench1.inp > mio1_15_s_16_1_4_4_d.dat

       $mycommand -np 16 ./main_2c_15_d.x 2 2 4 bench.inp > mio2_15_s_16_2_2_4_d.dat
       $mycommand -np 16 ./main_2c_15_d.x 1 4 4 bench.inp > mio2_15_s_16_1_4_4_d.dat

fi

if [[ "$mynodes" == "8" ]]; then

    echo "Mode: WEAK - Running weak configuration"
    
    $mycommand -np 32 ./main_1c_15_d.x 1 1 32 bench1w.inp > mio1_15_w_32_1_1_32_d.dat
    $mycommand -np 32 ./main_1c_15_d.x 1 4 8 bench1w.inp > mio1_15_w_32_1_4_8_d.dat
    $mycommand -np 32 ./main_1c_15_d.x 2 4 4 bench1w.inp > mio1_15_w_32_2_4_4_d.dat
    
    $mycommand -np 32 ./main_2c_15_d.x 1 1 32 benchw.inp > mio2_15_w_32_1_1_32_d.dat
    $mycommand -np 32 ./main_2c_15_d.x 1 4 8 benchw.inp > mio2_15_w_32_1_4_8_d.dat
    $mycommand -np 32 ./main_2c_15_d.x 2 4 4 benchw.inp > mio2_15_w_32_2_4_4_d.dat
    
    echo "Mode: STRONG - Running strong configuration"
    $mycommand -np 32 ./main_1c_15_d.x 2 4 4 bench1.inp > mio1_15_s_32_2_4_4_d.dat
    $mycommand -np 32 ./main_1c_15_d.x 1 4 8 bench1.inp > mio1_15_s_32_1_4_8_d.dat
    
    $mycommand -np 32 ./main_2c_15_d.x 2 4 4 bench.inp > mio2_15_s_32_2_4_4_d.dat
    $mycommand -np 32 ./main_2c_15_d.x 1 4 8 bench.inp > mio2_15_s_32_1_4_8_d.dat

fi

if [[ "$mynodes" == "16" ]]; then
    echo "Mode: WEAK - Running weak configuration"
    $mycommand -np 64 ./main_1c_15_d.x  1 1 64 bench1w.inp > mio1_15_w_64_1_1_64_d.dat
    $mycommand -np 64 ./main_1c_15_d.x  1 8 8 bench1w.inp > mio1_15_w_64_1_8_8_d.dat
    $mycommand -np 64 ./main_1c_15_d.x  2 4 8 bench1w.inp > mio1_15_w_64_2_4_8_d.dat
    $mycommand -np 64 ./main_1c_15_d.x  4 4 4 bench1w.inp > mio1_15_w_64_4_4_4_d.dat
    
    $mycommand -np 64 ./main_2c_15_d.x  1 1 64 benchw.inp > mio2_15_w_64_1_1_64_d.dat
    $mycommand -np 64 ./main_2c_15_d.x  1 8 8 benchw.inp > mio2_15_w_64_1_8_8_d.dat
    $mycommand -np 64 ./main_2c_15_d.x  2 4 8 benchw.inp > mio2_15_w_64_2_4_8_d.dat
    $mycommand -np 64 ./main_2c_15_d.x  4 4 4 benchw.inp > mio2_15_w_64_4_4_4_d.dat
    
    echo "Mode: STRONG - Running strong configuration"
    $mycommand -np 64 ./main_1c_15_d.x 4 4 4 bench1.inp > mio1_15_s_64_4_4_4_d.dat
    $mycommand -np 64 ./main_1c_15_d.x 2 4 8 bench1.inp > mio1_15_s_64_2_4_8_d.dat
    $mycommand -np 64 ./main_1c_15_d.x 1 8 8 bench1.inp > mio1_15_s_64_1_8_8_d.dat
    
    $mycommand -np 64 ./main_2c_15_d.x 4 4 4 bench.inp > mio2_15_s_64_4_4_4_d.dat
    $mycommand -np 64 ./main_2c_15_d.x 2 4 8 bench.inp > mio2_15_s_64_2_4_8_d.dat
    $mycommand -np 64 ./main_2c_15_d.x 1 8 8 bench.inp > mio2_15_s_64_1_8_8_d.dat

fi

if [[ "$mynodes" == "1" ]]; then
    echo "Mode: WEAK - Running weak configuration"
    $mycommand -np 1 ./main_1c_19_d.x 1 1 1 bench1w.inp > mio1_19_w_1_1_1_1_d.dat
    $mycommand -np 2 ./main_1c_19_d.x 1 1 2 bench1w.inp > mio1_19_w_2_1_1_2_d.dat
    $mycommand -np 4 ./main_1c_19_d.x 1 2 2 bench1w.inp > mio1_19_w_4_1_2_2_d.dat
    $mycommand -np 4 ./main_1c_19_d.x 1 1 4 bench1w.inp > mio1_19_w_4_1_1_4_d.dat

    $mycommand -np 1 ./main_2c_19_d.x 1 1 1 benchw.inp > mio2_19_w_1_1_1_1_d.dat
    $mycommand -np 2 ./main_2c_19_d.x 1 1 2 benchw.inp > mio2_19_w_2_1_1_2_d.dat
    $mycommand -np 4 ./main_2c_19_d.x 1 2 2 benchw.inp > mio2_19_w_4_1_2_2_d.dat
    $mycommand -np 4 ./main_2c_19_d.x 1 1 4 benchw.inp > mio2_19_w_4_1_1_4_d.dat
    
    echo "Mode: STRONG - Running strong configuration"
    
    $mycommand  -np 1 ./main_1c_19_d.x 1 1 1 bench1.inp > mio1_19_s_1_1_1_1_d.dat
    $mycommand  -np 2 ./main_1c_19_d.x 1 1 2 bench1.inp > mio1_19_s_2_1_1_2_d.dat
    $mycommand  -np 4 ./main_1c_19_d.x 1 1 4 bench1.inp > mio1_19_s_4_1_1_4_d.dat
    $mycommand  -np 4 ./main_1c_19_d.x 1 2 2 bench1.inp > mio1_19_s_4_1_2_2_d.dat

    $mycommand  -np 1 ./main_2c_19_d.x 1 1 1 bench.inp > mio2_19_s_1_1_1_1_d.dat
    $mycommand  -np 2 ./main_2c_19_d.x 1 1 2 bench.inp > mio2_19_s_2_1_1_2_d.dat
    $mycommand  -np 4 ./main_2c_19_d.x 1 1 4 bench.inp > mio2_19_s_4_1_1_4_d.dat
    $mycommand  -np 4 ./main_2c_19_d.x 1 2 2 bench.inp > mio2_19_s_4_1_2_2_d.dat
fi


if [[ "$mynodes" == "2" ]]; then
	echo "Mode: WEAK - Running weak configuration"
        $mycommand -np 8 ./main_1c_19_d.x 1 1 8 bench1w.inp > mio1_19_w_8_1_1_8_d.dat
        $mycommand -np 8 ./main_1c_19_d.x 1 2 4 bench1w.inp > mio1_19_w_8_1_2_4_d.dat
        $mycommand -np 8 ./main_1c_19_d.x 2 2 2 bench1w.inp > mio1_19_w_8_2_2_2_d.dat

	    $mycommand -np 8 ./main_2c_19_d.x 1 1 8 benchw.inp > mio2_19_w_8_1_1_8_d.dat
	    $mycommand -np 8 ./main_2c_19_d.x 1 2 4 benchw.inp > mio2_19_w_8_1_2_4_d.dat
        $mycommand -np 8 ./main_2c_19_d.x 2 2 2 benchw.inp > mio2_19_w_8_2_2_2_d.dat
	    
        echo "Mode: STRONG - Running strong configuration"
        $mycommand -np 8 ./main_1c_19_d.x  2 2 2 bench1.inp > mio1_19_s_8_2_2_2_d.dat
        $mycommand -np 8 ./main_1c_19_d.x  1 2 4 bench1.inp > mio1_19_s_8_1_2_4_d.dat

        $mycommand -np 8 ./main_2c_19_d.x  2 2 2 bench.inp > mio2_19_s_8_2_2_2_d.dat
        $mycommand -np 8 ./main_2c_19_d.x  1 2 4 bench.inp > mio2_19_s_8_1_2_4_d.dat

fi

if [[ "$mynodes" == "4" ]]; then
        echo "Mode: WEAK - Running weak configuration"

        $mycommand -np 16 ./main_1c_19_d.x 1 1 16 bench1w.inp > mio1_19_w_16_1_1_16_d.dat
        $mycommand -np 16 ./main_1c_19_d.x 1 4 4 bench1w.inp > mio1_19_w_16_1_4_4_d.dat
        $mycommand -np 16 ./main_1c_19_d.x 2 2 4 bench1w.inp > mio1_19_w_16_2_2_4_d.dat

	    $mycommand -np 16 ./main_2c_19_d.x 1 1 16 benchw.inp > mio2_19_w_16_1_1_16_d.dat
	    $mycommand -np 16 ./main_2c_19_d.x 1 4 4 benchw.inp > mio2_19_w_16_1_4_4_d.dat
        $mycommand -np 16 ./main_2c_19_d.x 2 2 4 benchw.inp > mio2_19_w_16_2_2_4_d.dat
	    
        echo "Mode: STRONG - Running strong configuration"
       $mycommand -np 16 ./main_1c_19_d.x 2 2 4 bench1.inp > mio1_19_s_16_2_2_4_d.dat
       $mycommand -np 16 ./main_1c_19_d.x 1 4 4 bench1.inp > mio1_19_s_16_1_4_4_d.dat

       $mycommand -np 16 ./main_2c_19_d.x 2 2 4 bench.inp > mio2_19_s_16_2_2_4_d.dat
       $mycommand -np 16 ./main_2c_19_d.x 1 4 4 bench.inp > mio2_19_s_16_1_4_4_d.dat

fi

if [[ "$mynodes" == "8" ]]; then

    echo "Mode: WEAK - Running weak configuration"
    
    $mycommand -np 32 ./main_1c_19_d.x 1 1 32 bench1w.inp > mio1_19_w_32_1_1_32_d.dat
    $mycommand -np 32 ./main_1c_19_d.x 1 4 8 bench1w.inp > mio1_19_w_32_1_4_8_d.dat
    $mycommand -np 32 ./main_1c_19_d.x 2 4 4 bench1w.inp > mio1_19_w_32_2_4_4_d.dat
    
    $mycommand -np 32 ./main_2c_19_d.x 1 1 32 benchw.inp > mio2_19_w_32_1_1_32_d.dat
    $mycommand -np 32 ./main_2c_19_d.x 1 4 8 benchw.inp > mio2_19_w_32_1_4_8_d.dat
    $mycommand -np 32 ./main_2c_19_d.x 2 4 4 benchw.inp > mio2_19_w_32_2_4_4_d.dat
    
    echo "Mode: STRONG - Running strong configuration"
    $mycommand -np 32 ./main_1c_19_d.x 2 4 4 bench1.inp > mio1_19_s_32_2_4_4_d.dat
    $mycommand -np 32 ./main_1c_19_d.x 1 4 8 bench1.inp > mio1_19_s_32_1_4_8_d.dat
    
    $mycommand -np 32 ./main_2c_19_d.x 2 4 4 bench.inp > mio2_19_s_32_2_4_4_d.dat
    $mycommand -np 32 ./main_2c_19_d.x 1 4 8 bench.inp > mio2_19_s_32_1_4_8_d.dat

fi

if [[ "$mynodes" == "16" ]]; then
    echo "Mode: WEAK - Running weak configuration"
    $mycommand -np 64 ./main_1c_19_d.x  1 1 64 bench1w.inp > mio1_19_w_64_1_1_64_d.dat
    $mycommand -np 64 ./main_1c_19_d.x  1 8 8 bench1w.inp > mio1_19_w_64_1_8_8_d.dat
    $mycommand -np 64 ./main_1c_19_d.x  2 4 8 bench1w.inp > mio1_19_w_64_2_4_8_d.dat
    $mycommand -np 64 ./main_1c_19_d.x  4 4 4 bench1w.inp > mio1_19_w_64_4_4_4_d.dat
    
    $mycommand -np 64 ./main_2c_19_d.x  1 1 64 benchw.inp > mio2_19_w_64_1_1_64_d.dat
    $mycommand -np 64 ./main_2c_19_d.x  1 8 8 benchw.inp > mio2_19_w_64_1_8_8_d.dat
    $mycommand -np 64 ./main_2c_19_d.x  2 4 8 benchw.inp > mio2_19_w_64_2_4_8_d.dat
    $mycommand -np 64 ./main_2c_19_d.x  4 4 4 benchw.inp > mio2_19_w_64_4_4_4_d.dat
    
    echo "Mode: STRONG - Running strong configuration"
    $mycommand -np 64 ./main_1c_19_d.x 4 4 4 bench1.inp > mio1_19_s_64_4_4_4_d.dat
    $mycommand -np 64 ./main_1c_19_d.x 2 4 8 bench1.inp > mio1_19_s_64_2_4_8_d.dat
    $mycommand -np 64 ./main_1c_19_d.x 1 8 8 bench1.inp > mio1_19_s_64_1_8_8_d.dat
    
    $mycommand -np 64 ./main_2c_19_d.x 4 4 4 bench.inp > mio2_19_s_64_4_4_4_d.dat
    $mycommand -np 64 ./main_2c_19_d.x 2 4 8 bench.inp > mio2_19_s_64_2_4_8_d.dat
    $mycommand -np 64 ./main_2c_19_d.x 1 8 8 bench.inp > mio2_19_s_64_1_8_8_d.dat

fi


if [[ "$mynodes" == "1" ]]; then
    echo "Mode: WEAK - Running weak configuration"
    $mycommand -np 1 ./main_1c_27_d.x 1 1 1 bench1w.inp > mio1_27_w_1_1_1_1_d.dat
    $mycommand -np 2 ./main_1c_27_d.x 1 1 2 bench1w.inp > mio1_27_w_2_1_1_2_d.dat
    $mycommand -np 4 ./main_1c_27_d.x 1 2 2 bench1w.inp > mio1_27_w_4_1_2_2_d.dat
    $mycommand -np 4 ./main_1c_27_d.x 1 1 4 bench1w.inp > mio1_27_w_4_1_1_4_d.dat

    $mycommand -np 1 ./main_2c_27_d.x 1 1 1 benchw.inp > mio2_27_w_1_1_1_1_d.dat
    $mycommand -np 2 ./main_2c_27_d.x 1 1 2 benchw.inp > mio2_27_w_2_1_1_2_d.dat
    $mycommand -np 4 ./main_2c_27_d.x 1 2 2 benchw.inp > mio2_27_w_4_1_2_2_d.dat
    $mycommand -np 4 ./main_2c_27_d.x 1 1 4 benchw.inp > mio2_27_w_4_1_1_4_d.dat
    
    echo "Mode: STRONG - Running strong configuration"
    
    $mycommand  -np 1 ./main_1c_27_d.x 1 1 1 bench1.inp > mio1_27_s_1_1_1_1_d.dat
    $mycommand  -np 2 ./main_1c_27_d.x 1 1 2 bench1.inp > mio1_27_s_2_1_1_2_d.dat
    $mycommand  -np 4 ./main_1c_27_d.x 1 1 4 bench1.inp > mio1_27_s_4_1_1_4_d.dat
    $mycommand  -np 4 ./main_1c_27_d.x 1 2 2 bench1.inp > mio1_27_s_4_1_2_2_d.dat

    $mycommand  -np 1 ./main_2c_27_d.x 1 1 1 bench.inp > mio2_27_s_1_1_1_1_d.dat
    $mycommand  -np 2 ./main_2c_27_d.x 1 1 2 bench.inp > mio2_27_s_2_1_1_2_d.dat
    $mycommand  -np 4 ./main_2c_27_d.x 1 1 4 bench.inp > mio2_27_s_4_1_1_4_d.dat
    $mycommand  -np 4 ./main_2c_27_d.x 1 2 2 bench.inp > mio2_27_s_4_1_2_2_d.dat
fi


if [[ "$mynodes" == "2" ]]; then
	echo "Mode: WEAK - Running weak configuration"
        $mycommand -np 8 ./main_1c_27_d.x 1 1 8 bench1w.inp > mio1_27_w_8_1_1_8_d.dat
        $mycommand -np 8 ./main_1c_27_d.x 1 2 4 bench1w.inp > mio1_27_w_8_1_2_4_d.dat
        $mycommand -np 8 ./main_1c_27_d.x 2 2 2 bench1w.inp > mio1_27_w_8_2_2_2_d.dat

	    $mycommand -np 8 ./main_2c_27_d.x 1 1 8 benchw.inp > mio2_27_w_8_1_1_8_d.dat
	    $mycommand -np 8 ./main_2c_27_d.x 1 2 4 benchw.inp > mio2_27_w_8_1_2_4_d.dat
        $mycommand -np 8 ./main_2c_27_d.x 2 2 2 benchw.inp > mio2_27_w_8_2_2_2_d.dat
	    
        echo "Mode: STRONG - Running strong configuration"
        $mycommand -np 8 ./main_1c_27_d.x  2 2 2 bench1.inp > mio1_27_s_8_2_2_2_d.dat
        $mycommand -np 8 ./main_1c_27_d.x  1 2 4 bench1.inp > mio1_27_s_8_1_2_4_d.dat

        $mycommand -np 8 ./main_2c_27_d.x  2 2 2 bench.inp > mio2_27_s_8_2_2_2_d.dat
        $mycommand -np 8 ./main_2c_27_d.x  1 2 4 bench.inp > mio2_27_s_8_1_2_4_d.dat

fi

if [[ "$mynodes" == "4" ]]; then
        echo "Mode: WEAK - Running weak configuration"

        $mycommand -np 16 ./main_1c_27_d.x 1 1 16 bench1w.inp > mio1_27_w_16_1_1_16_d.dat
        $mycommand -np 16 ./main_1c_27_d.x 1 4 4 bench1w.inp > mio1_27_w_16_1_4_4_d.dat
        $mycommand -np 16 ./main_1c_27_d.x 2 2 4 bench1w.inp > mio1_27_w_16_2_2_4_d.dat

	    $mycommand -np 16 ./main_2c_27_d.x 1 1 16 benchw.inp > mio2_27_w_16_1_1_16_d.dat
	    $mycommand -np 16 ./main_2c_27_d.x 1 4 4 benchw.inp > mio2_27_w_16_1_4_4_d.dat
        $mycommand -np 16 ./main_2c_27_d.x 2 2 4 benchw.inp > mio2_27_w_16_2_2_4_d.dat
	    
        echo "Mode: STRONG - Running strong configuration"
       $mycommand -np 16 ./main_1c_27_d.x 2 2 4 bench1.inp > mio1_27_s_16_2_2_4_d.dat
       $mycommand -np 16 ./main_1c_27_d.x 1 4 4 bench1.inp > mio1_27_s_16_1_4_4_d.dat

       $mycommand -np 16 ./main_2c_27_d.x 2 2 4 bench.inp > mio2_27_s_16_2_2_4_d.dat
       $mycommand -np 16 ./main_2c_27_d.x 1 4 4 bench.inp > mio2_27_s_16_1_4_4_d.dat

fi

if [[ "$mynodes" == "8" ]]; then

    echo "Mode: WEAK - Running weak configuration"
    
    $mycommand -np 32 ./main_1c_27_d.x 1 1 32 bench1w.inp > mio1_27_w_32_1_1_32_d.dat
    $mycommand -np 32 ./main_1c_27_d.x 1 4 8 bench1w.inp > mio1_27_w_32_1_4_8_d.dat
    $mycommand -np 32 ./main_1c_27_d.x 2 4 4 bench1w.inp > mio1_27_w_32_2_4_4_d.dat
    
    $mycommand -np 32 ./main_2c_27_d.x 1 1 32 benchw.inp > mio2_27_w_32_1_1_32_d.dat
    $mycommand -np 32 ./main_2c_27_d.x 1 4 8 benchw.inp > mio2_27_w_32_1_4_8_d.dat
    $mycommand -np 32 ./main_2c_27_d.x 2 4 4 benchw.inp > mio2_27_w_32_2_4_4_d.dat
    
    echo "Mode: STRONG - Running strong configuration"
    $mycommand -np 32 ./main_1c_27_d.x 2 4 4 bench1.inp > mio1_27_s_32_2_4_4_d.dat
    $mycommand -np 32 ./main_1c_27_d.x 1 4 8 bench1.inp > mio1_27_s_32_1_4_8_d.dat
    
    $mycommand -np 32 ./main_2c_27_d.x 2 4 4 bench.inp > mio2_27_s_32_2_4_4_d.dat
    $mycommand -np 32 ./main_2c_27_d.x 1 4 8 bench.inp > mio2_27_s_32_1_4_8_d.dat

fi

if [[ "$mynodes" == "16" ]]; then
    echo "Mode: WEAK - Running weak configuration"
    $mycommand -np 64 ./main_1c_27_d.x  1 1 64 bench1w.inp > mio1_27_w_64_1_1_64_d.dat
    $mycommand -np 64 ./main_1c_27_d.x  1 8 8 bench1w.inp > mio1_27_w_64_1_8_8_d.dat
    $mycommand -np 64 ./main_1c_27_d.x  2 4 8 bench1w.inp > mio1_27_w_64_2_4_8_d.dat
    $mycommand -np 64 ./main_1c_27_d.x  4 4 4 bench1w.inp > mio1_27_w_64_4_4_4_d.dat
    
    $mycommand -np 64 ./main_2c_27_d.x  1 1 64 benchw.inp > mio2_27_w_64_1_1_64_d.dat
    $mycommand -np 64 ./main_2c_27_d.x  1 8 8 benchw.inp > mio2_27_w_64_1_8_8_d.dat
    $mycommand -np 64 ./main_2c_27_d.x  2 4 8 benchw.inp > mio2_27_w_64_2_4_8_d.dat
    $mycommand -np 64 ./main_2c_27_d.x  4 4 4 benchw.inp > mio2_27_w_64_4_4_4_d.dat
    
    echo "Mode: STRONG - Running strong configuration"
    $mycommand -np 64 ./main_1c_27_d.x 4 4 4 bench1.inp > mio1_27_s_64_4_4_4_d.dat
    $mycommand -np 64 ./main_1c_27_d.x 2 4 8 bench1.inp > mio1_27_s_64_2_4_8_d.dat
    $mycommand -np 64 ./main_1c_27_d.x 1 8 8 bench1.inp > mio1_27_s_64_1_8_8_d.dat
    
    $mycommand -np 64 ./main_2c_27_d.x 4 4 4 bench.inp > mio2_27_s_64_4_4_4_d.dat
    $mycommand -np 64 ./main_2c_27_d.x 2 4 8 bench.inp > mio2_27_s_64_2_4_8_d.dat
    $mycommand -np 64 ./main_2c_27_d.x 1 8 8 bench.inp > mio2_27_s_64_1_8_8_d.dat

fi


if [[ "$mynodes" == "1" ]]; then
    echo "Mode: WEAK - Running weak configuration"
    $mycommand -np 1 ./main_1c_27high_d.x 1 1 1 bench1w.inp > mio1_27h_w_1_1_1_1_d.dat
    $mycommand -np 2 ./main_1c_27high_d.x 1 1 2 bench1w.inp > mio1_27h_w_2_1_1_2_d.dat
    $mycommand -np 4 ./main_1c_27high_d.x 1 2 2 bench1w.inp > mio1_27h_w_4_1_2_2_d.dat
    $mycommand -np 4 ./main_1c_27high_d.x 1 1 4 bench1w.inp > mio1_27h_w_4_1_1_4_d.dat

    $mycommand -np 1 ./main_2c_27high_d.x 1 1 1 benchw.inp > mio2_27h_w_1_1_1_1_d.dat
    $mycommand -np 2 ./main_2c_27high_d.x 1 1 2 benchw.inp > mio2_27h_w_2_1_1_2_d.dat
    $mycommand -np 4 ./main_2c_27high_d.x 1 2 2 benchw.inp > mio2_27h_w_4_1_2_2_d.dat
    $mycommand -np 4 ./main_2c_27high_d.x 1 1 4 benchw.inp > mio2_27h_w_4_1_1_4_d.dat
    
    echo "Mode: STRONG - Running strong configuration"
    
    $mycommand  -np 1 ./main_1c_27high_d.x 1 1 1 bench1.inp > mio1_27h_s_1_1_1_1_d.dat
    $mycommand  -np 2 ./main_1c_27high_d.x 1 1 2 bench1.inp > mio1_27h_s_2_1_1_2_d.dat
    $mycommand  -np 4 ./main_1c_27high_d.x 1 1 4 bench1.inp > mio1_27h_s_4_1_1_4_d.dat
    $mycommand  -np 4 ./main_1c_27high_d.x 1 2 2 bench1.inp > mio1_27h_s_4_1_2_2_d.dat

    $mycommand  -np 1 ./main_2c_27high_d.x 1 1 1 bench.inp > mio2_27h_s_1_1_1_1_d.dat
    $mycommand  -np 2 ./main_2c_27high_d.x 1 1 2 bench.inp > mio2_27h_s_2_1_1_2_d.dat
    $mycommand  -np 4 ./main_2c_27high_d.x 1 1 4 bench.inp > mio2_27h_s_4_1_1_4_d.dat
    $mycommand  -np 4 ./main_2c_27high_d.x 1 2 2 bench.inp > mio2_27h_s_4_1_2_2_d.dat
fi


if [[ "$mynodes" == "2" ]]; then
	echo "Mode: WEAK - Running weak configuration"
        $mycommand -np 8 ./main_1c_27high_d.x 1 1 8 bench1w.inp > mio1_27h_w_8_1_1_8_d.dat
        $mycommand -np 8 ./main_1c_27high_d.x 1 2 4 bench1w.inp > mio1_27h_w_8_1_2_4_d.dat
        $mycommand -np 8 ./main_1c_27high_d.x 2 2 2 bench1w.inp > mio1_27h_w_8_2_2_2_d.dat

	    $mycommand -np 8 ./main_2c_27high_d.x 1 1 8 benchw.inp > mio2_27h_w_8_1_1_8_d.dat
	    $mycommand -np 8 ./main_2c_27high_d.x 1 2 4 benchw.inp > mio2_27h_w_8_1_2_4_d.dat
        $mycommand -np 8 ./main_2c_27high_d.x 2 2 2 benchw.inp > mio2_27h_w_8_2_2_2_d.dat
	    
        echo "Mode: STRONG - Running strong configuration"
        $mycommand -np 8 ./main_1c_27high_d.x  2 2 2 bench1.inp > mio1_27h_s_8_2_2_2_d.dat
        $mycommand -np 8 ./main_1c_27high_d.x  1 2 4 bench1.inp > mio1_27h_s_8_1_2_4_d.dat

        $mycommand -np 8 ./main_2c_27high_d.x  2 2 2 bench.inp > mio2_27h_s_8_2_2_2_d.dat
        $mycommand -np 8 ./main_2c_27high_d.x  1 2 4 bench.inp > mio2_27h_s_8_1_2_4_d.dat

fi

if [[ "$mynodes" == "4" ]]; then
        echo "Mode: WEAK - Running weak configuration"

        $mycommand -np 16 ./main_1c_27high_d.x 1 1 16 bench1w.inp > mio1_27h_w_16_1_1_16_d.dat
        $mycommand -np 16 ./main_1c_27high_d.x 1 4 4 bench1w.inp > mio1_27h_w_16_1_4_4_d.dat
        $mycommand -np 16 ./main_1c_27high_d.x 2 2 4 bench1w.inp > mio1_27h_w_16_2_2_4_d.dat

	    $mycommand -np 16 ./main_2c_27high_d.x 1 1 16 benchw.inp > mio2_27h_w_16_1_1_16_d.dat
	    $mycommand -np 16 ./main_2c_27high_d.x 1 4 4 benchw.inp > mio2_27h_w_16_1_4_4_d.dat
        $mycommand -np 16 ./main_2c_27high_d.x 2 2 4 benchw.inp > mio2_27h_w_16_2_2_4_d.dat
	    
        echo "Mode: STRONG - Running strong configuration"
       $mycommand -np 16 ./main_1c_27high_d.x 2 2 4 bench1.inp > mio1_27h_s_16_2_2_4_d.dat
       $mycommand -np 16 ./main_1c_27high_d.x 1 4 4 bench1.inp > mio1_27h_s_16_1_4_4_d.dat

       $mycommand -np 16 ./main_2c_27high_d.x 2 2 4 bench.inp > mio2_27h_s_16_2_2_4_d.dat
       $mycommand -np 16 ./main_2c_27high_d.x 1 4 4 bench.inp > mio2_27h_s_16_1_4_4_d.dat

fi

if [[ "$mynodes" == "8" ]]; then

    echo "Mode: WEAK - Running weak configuration"
    
    $mycommand -np 32 ./main_1c_27high_d.x 1 1 32 bench1w.inp > mio1_27h_w_32_1_1_32_d.dat
    $mycommand -np 32 ./main_1c_27high_d.x 1 4 8 bench1w.inp > mio1_27h_w_32_1_4_8_d.dat
    $mycommand -np 32 ./main_1c_27high_d.x 2 4 4 bench1w.inp > mio1_27h_w_32_2_4_4_d.dat
    
    $mycommand -np 32 ./main_2c_27high_d.x 1 1 32 benchw.inp > mio2_27h_w_32_1_1_32_d.dat
    $mycommand -np 32 ./main_2c_27high_d.x 1 4 8 benchw.inp > mio2_27h_w_32_1_4_8_d.dat
    $mycommand -np 32 ./main_2c_27high_d.x 2 4 4 benchw.inp > mio2_27h_w_32_2_4_4_d.dat
    
    echo "Mode: STRONG - Running strong configuration"
    $mycommand -np 32 ./main_1c_27high_d.x 2 4 4 bench1.inp > mio1_27h_s_32_2_4_4_d.dat
    $mycommand -np 32 ./main_1c_27high_d.x 1 4 8 bench1.inp > mio1_27h_s_32_1_4_8_d.dat
    
    $mycommand -np 32 ./main_2c_27high_d.x 2 4 4 bench.inp > mio2_27h_s_32_2_4_4_d.dat
    $mycommand -np 32 ./main_2c_27high_d.x 1 4 8 bench.inp > mio2_27h_s_32_1_4_8_d.dat

fi

if [[ "$mynodes" == "16" ]]; then
    echo "Mode: WEAK - Running weak configuration"
    $mycommand -np 64 ./main_1c_27high_d.x  1 1 64 bench1w.inp > mio1_27h_w_64_1_1_64_d.dat
    $mycommand -np 64 ./main_1c_27high_d.x  1 8 8 bench1w.inp > mio1_27h_w_64_1_8_8_d.dat
    $mycommand -np 64 ./main_1c_27high_d.x  2 4 8 bench1w.inp > mio1_27h_w_64_2_4_8_d.dat
    $mycommand -np 64 ./main_1c_27high_d.x  4 4 4 bench1w.inp > mio1_27h_w_64_4_4_4_d.dat
    
    $mycommand -np 64 ./main_2c_27high_d.x  1 1 64 benchw.inp > mio2_27h_w_64_1_1_64_d.dat
    $mycommand -np 64 ./main_2c_27high_d.x  1 8 8 benchw.inp > mio2_27h_w_64_1_8_8_d.dat
    $mycommand -np 64 ./main_2c_27high_d.x  2 4 8 benchw.inp > mio2_27h_w_64_2_4_8_d.dat
    $mycommand -np 64 ./main_2c_27high_d.x  4 4 4 benchw.inp > mio2_27h_w_64_4_4_4_d.dat
    
    echo "Mode: STRONG - Running strong configuration"
    $mycommand -np 64 ./main_1c_27high_d.x 4 4 4 bench1.inp > mio1_27h_s_64_4_4_4_d.dat
    $mycommand -np 64 ./main_1c_27high_d.x 2 4 8 bench1.inp > mio1_27h_s_64_2_4_8_d.dat
    $mycommand -np 64 ./main_1c_27high_d.x 1 8 8 bench1.inp > mio1_27h_s_64_1_8_8_d.dat
    
    $mycommand -np 64 ./main_2c_27high_d.x 4 4 4 bench.inp > mio2_27h_s_64_4_4_4_d.dat
    $mycommand -np 64 ./main_2c_27high_d.x 2 4 8 bench.inp > mio2_27h_s_64_2_4_8_d.dat
    $mycommand -np 64 ./main_2c_27high_d.x 1 8 8 bench.inp > mio2_27h_s_64_1_8_8_d.dat

fi


###########################


if [[ "$mynodes" == "1" ]]; then
    echo "Mode: WEAK - Running weak configuration"
    $mycommand -np 1 ./main_1c_15_sd.x 1 1 1 bench1w.inp > mio1_15_w_1_1_1_1_sd.dat
    $mycommand -np 2 ./main_1c_15_sd.x 1 1 2 bench1w.inp > mio1_15_w_2_1_1_2_sd.dat
    $mycommand -np 4 ./main_1c_15_sd.x 1 2 2 bench1w.inp > mio1_15_w_4_1_2_2_sd.dat
    $mycommand -np 4 ./main_1c_15_sd.x 1 1 4 bench1w.inp > mio1_15_w_4_1_1_4_sd.dat

    $mycommand -np 1 ./main_2c_15_sd.x 1 1 1 benchw.inp > mio2_15_w_1_1_1_1_sd.dat
    $mycommand -np 2 ./main_2c_15_sd.x 1 1 2 benchw.inp > mio2_15_w_2_1_1_2_sd.dat
    $mycommand -np 4 ./main_2c_15_sd.x 1 2 2 benchw.inp > mio2_15_w_4_1_2_2_sd.dat
    $mycommand -np 4 ./main_2c_15_sd.x 1 1 4 benchw.inp > mio2_15_w_4_1_1_4_sd.dat
    
    echo "Mode: STRONG - Running strong configuration"
    
    $mycommand  -np 1 ./main_1c_15_sd.x 1 1 1 bench1.inp > mio1_15_s_1_1_1_1_sd.dat
    $mycommand  -np 2 ./main_1c_15_sd.x 1 1 2 bench1.inp > mio1_15_s_2_1_1_2_sd.dat
    $mycommand  -np 4 ./main_1c_15_sd.x 1 1 4 bench1.inp > mio1_15_s_4_1_1_4_sd.dat
    $mycommand  -np 4 ./main_1c_15_sd.x 1 2 2 bench1.inp > mio1_15_s_4_1_2_2_sd.dat

    $mycommand  -np 1 ./main_2c_15_sd.x 1 1 1 bench.inp > mio2_15_s_1_1_1_1_sd.dat
    $mycommand  -np 2 ./main_2c_15_sd.x 1 1 2 bench.inp > mio2_15_s_2_1_1_2_sd.dat
    $mycommand  -np 4 ./main_2c_15_sd.x 1 1 4 bench.inp > mio2_15_s_4_1_1_4_sd.dat
    $mycommand  -np 4 ./main_2c_15_sd.x 1 2 2 bench.inp > mio2_15_s_4_1_2_2_sd.dat
fi


if [[ "$mynodes" == "2" ]]; then
	echo "Mode: WEAK - Running weak configuration"
        $mycommand -np 8 ./main_1c_15_sd.x 1 1 8 bench1w.inp > mio1_15_w_8_1_1_8_sd.dat
        $mycommand -np 8 ./main_1c_15_sd.x 1 2 4 bench1w.inp > mio1_15_w_8_1_2_4_sd.dat
        $mycommand -np 8 ./main_1c_15_sd.x 2 2 2 bench1w.inp > mio1_15_w_8_2_2_2_sd.dat

	    $mycommand -np 8 ./main_2c_15_sd.x 1 1 8 benchw.inp > mio2_15_w_8_1_1_8_sd.dat
	    $mycommand -np 8 ./main_2c_15_sd.x 1 2 4 benchw.inp > mio2_15_w_8_1_2_4_sd.dat
        $mycommand -np 8 ./main_2c_15_sd.x 2 2 2 benchw.inp > mio2_15_w_8_2_2_2_sd.dat
	    
        echo "Mode: STRONG - Running strong configuration"
        $mycommand -np 8 ./main_1c_15_sd.x  2 2 2 bench1.inp > mio1_15_s_8_2_2_2_sd.dat
        $mycommand -np 8 ./main_1c_15_sd.x  1 2 4 bench1.inp > mio1_15_s_8_1_2_4_sd.dat

        $mycommand -np 8 ./main_2c_15_sd.x  2 2 2 bench.inp > mio2_15_s_8_2_2_2_sd.dat
        $mycommand -np 8 ./main_2c_15_sd.x  1 2 4 bench.inp > mio2_15_s_8_1_2_4_sd.dat

fi

if [[ "$mynodes" == "4" ]]; then
        echo "Mode: WEAK - Running weak configuration"

        $mycommand -np 16 ./main_1c_15_sd.x 1 1 16 bench1w.inp > mio1_15_w_16_1_1_16_sd.dat
        $mycommand -np 16 ./main_1c_15_sd.x 1 4 4 bench1w.inp > mio1_15_w_16_1_4_4_sd.dat
        $mycommand -np 16 ./main_1c_15_sd.x 2 2 4 bench1w.inp > mio1_15_w_16_2_2_4_sd.dat

	    $mycommand -np 16 ./main_2c_15_sd.x 1 1 16 benchw.inp > mio2_15_w_16_1_1_16_sd.dat
	    $mycommand -np 16 ./main_2c_15_sd.x 1 4 4 benchw.inp > mio2_15_w_16_1_4_4_sd.dat
        $mycommand -np 16 ./main_2c_15_sd.x 2 2 4 benchw.inp > mio2_15_w_16_2_2_4_sd.dat
	    
        echo "Mode: STRONG - Running strong configuration"
       $mycommand -np 16 ./main_1c_15_sd.x 2 2 4 bench1.inp > mio1_15_s_16_2_2_4_sd.dat
       $mycommand -np 16 ./main_1c_15_sd.x 1 4 4 bench1.inp > mio1_15_s_16_1_4_4_sd.dat

       $mycommand -np 16 ./main_2c_15_sd.x 2 2 4 bench.inp > mio2_15_s_16_2_2_4_sd.dat
       $mycommand -np 16 ./main_2c_15_sd.x 1 4 4 bench.inp > mio2_15_s_16_1_4_4_sd.dat

fi

if [[ "$mynodes" == "8" ]]; then

    echo "Mode: WEAK - Running weak configuration"
    
    $mycommand -np 32 ./main_1c_15_sd.x 1 1 32 bench1w.inp > mio1_15_w_32_1_1_32_sd.dat
    $mycommand -np 32 ./main_1c_15_sd.x 1 4 8 bench1w.inp > mio1_15_w_32_1_4_8_sd.dat
    $mycommand -np 32 ./main_1c_15_sd.x 2 4 4 bench1w.inp > mio1_15_w_32_2_4_4_sd.dat
    
    $mycommand -np 32 ./main_2c_15_sd.x 1 1 32 benchw.inp > mio2_15_w_32_1_1_32_sd.dat
    $mycommand -np 32 ./main_2c_15_sd.x 1 4 8 benchw.inp > mio2_15_w_32_1_4_8_sd.dat
    $mycommand -np 32 ./main_2c_15_sd.x 2 4 4 benchw.inp > mio2_15_w_32_2_4_4_sd.dat
    
    echo "Mode: STRONG - Running strong configuration"
    $mycommand -np 32 ./main_1c_15_sd.x 2 4 4 bench1.inp > mio1_15_s_32_2_4_4_sd.dat
    $mycommand -np 32 ./main_1c_15_sd.x 1 4 8 bench1.inp > mio1_15_s_32_1_4_8_sd.dat
    
    $mycommand -np 32 ./main_2c_15_sd.x 2 4 4 bench.inp > mio2_15_s_32_2_4_4_sd.dat
    $mycommand -np 32 ./main_2c_15_sd.x 1 4 8 bench.inp > mio2_15_s_32_1_4_8_sd.dat

fi

if [[ "$mynodes" == "16" ]]; then
    echo "Mode: WEAK - Running weak configuration"
    $mycommand -np 64 ./main_1c_15_sd.x  1 1 64 bench1w.inp > mio1_15_w_64_1_1_64_sd.dat
    $mycommand -np 64 ./main_1c_15_sd.x  1 8 8 bench1w.inp > mio1_15_w_64_1_8_8_sd.dat
    $mycommand -np 64 ./main_1c_15_sd.x  2 4 8 bench1w.inp > mio1_15_w_64_2_4_8_sd.dat
    $mycommand -np 64 ./main_1c_15_sd.x  4 4 4 bench1w.inp > mio1_15_w_64_4_4_4_sd.dat
    
    $mycommand -np 64 ./main_2c_15_sd.x  1 1 64 benchw.inp > mio2_15_w_64_1_1_64_sd.dat
    $mycommand -np 64 ./main_2c_15_sd.x  1 8 8 benchw.inp > mio2_15_w_64_1_8_8_sd.dat
    $mycommand -np 64 ./main_2c_15_sd.x  2 4 8 benchw.inp > mio2_15_w_64_2_4_8_sd.dat
    $mycommand -np 64 ./main_2c_15_sd.x  4 4 4 benchw.inp > mio2_15_w_64_4_4_4_sd.dat
    
    echo "Mode: STRONG - Running strong configuration"
    $mycommand -np 64 ./main_1c_15_sd.x 4 4 4 bench1.inp > mio1_15_s_64_4_4_4_sd.dat
    $mycommand -np 64 ./main_1c_15_sd.x 2 4 8 bench1.inp > mio1_15_s_64_2_4_8_sd.dat
    $mycommand -np 64 ./main_1c_15_sd.x 1 8 8 bench1.inp > mio1_15_s_64_1_8_8_sd.dat
    
    $mycommand -np 64 ./main_2c_15_sd.x 4 4 4 bench.inp > mio2_15_s_64_4_4_4_sd.dat
    $mycommand -np 64 ./main_2c_15_sd.x 2 4 8 bench.inp > mio2_15_s_64_2_4_8_sd.dat
    $mycommand -np 64 ./main_2c_15_sd.x 1 8 8 bench.inp > mio2_15_s_64_1_8_8_sd.dat

fi

if [[ "$mynodes" == "1" ]]; then
    echo "Mode: WEAK - Running weak configuration"
    $mycommand -np 1 ./main_1c_19_sd.x 1 1 1 bench1w.inp > mio1_19_w_1_1_1_1_sd.dat
    $mycommand -np 2 ./main_1c_19_sd.x 1 1 2 bench1w.inp > mio1_19_w_2_1_1_2_sd.dat
    $mycommand -np 4 ./main_1c_19_sd.x 1 2 2 bench1w.inp > mio1_19_w_4_1_2_2_sd.dat
    $mycommand -np 4 ./main_1c_19_sd.x 1 1 4 bench1w.inp > mio1_19_w_4_1_1_4_sd.dat

    $mycommand -np 1 ./main_2c_19_sd.x 1 1 1 benchw.inp > mio2_19_w_1_1_1_1_sd.dat
    $mycommand -np 2 ./main_2c_19_sd.x 1 1 2 benchw.inp > mio2_19_w_2_1_1_2_sd.dat
    $mycommand -np 4 ./main_2c_19_sd.x 1 2 2 benchw.inp > mio2_19_w_4_1_2_2_sd.dat
    $mycommand -np 4 ./main_2c_19_sd.x 1 1 4 benchw.inp > mio2_19_w_4_1_1_4_sd.dat
    
    echo "Mode: STRONG - Running strong configuration"
    
    $mycommand  -np 1 ./main_1c_19_sd.x 1 1 1 bench1.inp > mio1_19_s_1_1_1_1_sd.dat
    $mycommand  -np 2 ./main_1c_19_sd.x 1 1 2 bench1.inp > mio1_19_s_2_1_1_2_sd.dat
    $mycommand  -np 4 ./main_1c_19_sd.x 1 1 4 bench1.inp > mio1_19_s_4_1_1_4_sd.dat
    $mycommand  -np 4 ./main_1c_19_sd.x 1 2 2 bench1.inp > mio1_19_s_4_1_2_2_sd.dat

    $mycommand  -np 1 ./main_2c_19_sd.x 1 1 1 bench.inp > mio2_19_s_1_1_1_1_sd.dat
    $mycommand  -np 2 ./main_2c_19_sd.x 1 1 2 bench.inp > mio2_19_s_2_1_1_2_sd.dat
    $mycommand  -np 4 ./main_2c_19_sd.x 1 1 4 bench.inp > mio2_19_s_4_1_1_4_sd.dat
    $mycommand  -np 4 ./main_2c_19_sd.x 1 2 2 bench.inp > mio2_19_s_4_1_2_2_sd.dat
fi


if [[ "$mynodes" == "2" ]]; then
	echo "Mode: WEAK - Running weak configuration"
        $mycommand -np 8 ./main_1c_19_sd.x 1 1 8 bench1w.inp > mio1_19_w_8_1_1_8_sd.dat
        $mycommand -np 8 ./main_1c_19_sd.x 1 2 4 bench1w.inp > mio1_19_w_8_1_2_4_sd.dat
        $mycommand -np 8 ./main_1c_19_sd.x 2 2 2 bench1w.inp > mio1_19_w_8_2_2_2_sd.dat

	    $mycommand -np 8 ./main_2c_19_sd.x 1 1 8 benchw.inp > mio2_19_w_8_1_1_8_sd.dat
	    $mycommand -np 8 ./main_2c_19_sd.x 1 2 4 benchw.inp > mio2_19_w_8_1_2_4_sd.dat
        $mycommand -np 8 ./main_2c_19_sd.x 2 2 2 benchw.inp > mio2_19_w_8_2_2_2_sd.dat
	    
        echo "Mode: STRONG - Running strong configuration"
        $mycommand -np 8 ./main_1c_19_sd.x  2 2 2 bench1.inp > mio1_19_s_8_2_2_2_sd.dat
        $mycommand -np 8 ./main_1c_19_sd.x  1 2 4 bench1.inp > mio1_19_s_8_1_2_4_sd.dat

        $mycommand -np 8 ./main_2c_19_sd.x  2 2 2 bench.inp > mio2_19_s_8_2_2_2_sd.dat
        $mycommand -np 8 ./main_2c_19_sd.x  1 2 4 bench.inp > mio2_19_s_8_1_2_4_sd.dat

fi

if [[ "$mynodes" == "4" ]]; then
        echo "Mode: WEAK - Running weak configuration"

        $mycommand -np 16 ./main_1c_19_sd.x 1 1 16 bench1w.inp > mio1_19_w_16_1_1_16_sd.dat
        $mycommand -np 16 ./main_1c_19_sd.x 1 4 4 bench1w.inp > mio1_19_w_16_1_4_4_sd.dat
        $mycommand -np 16 ./main_1c_19_sd.x 2 2 4 bench1w.inp > mio1_19_w_16_2_2_4_sd.dat

	    $mycommand -np 16 ./main_2c_19_sd.x 1 1 16 benchw.inp > mio2_19_w_16_1_1_16_sd.dat
	    $mycommand -np 16 ./main_2c_19_sd.x 1 4 4 benchw.inp > mio2_19_w_16_1_4_4_sd.dat
        $mycommand -np 16 ./main_2c_19_sd.x 2 2 4 benchw.inp > mio2_19_w_16_2_2_4_sd.dat
	    
        echo "Mode: STRONG - Running strong configuration"
       $mycommand -np 16 ./main_1c_19_sd.x 2 2 4 bench1.inp > mio1_19_s_16_2_2_4_sd.dat
       $mycommand -np 16 ./main_1c_19_sd.x 1 4 4 bench1.inp > mio1_19_s_16_1_4_4_sd.dat

       $mycommand -np 16 ./main_2c_19_sd.x 2 2 4 bench.inp > mio2_19_s_16_2_2_4_sd.dat
       $mycommand -np 16 ./main_2c_19_sd.x 1 4 4 bench.inp > mio2_19_s_16_1_4_4_sd.dat

fi

if [[ "$mynodes" == "8" ]]; then

    echo "Mode: WEAK - Running weak configuration"
    
    $mycommand -np 32 ./main_1c_19_sd.x 1 1 32 bench1w.inp > mio1_19_w_32_1_1_32_sd.dat
    $mycommand -np 32 ./main_1c_19_sd.x 1 4 8 bench1w.inp > mio1_19_w_32_1_4_8_sd.dat
    $mycommand -np 32 ./main_1c_19_sd.x 2 4 4 bench1w.inp > mio1_19_w_32_2_4_4_sd.dat
    
    $mycommand -np 32 ./main_2c_19_sd.x 1 1 32 benchw.inp > mio2_19_w_32_1_1_32_sd.dat
    $mycommand -np 32 ./main_2c_19_sd.x 1 4 8 benchw.inp > mio2_19_w_32_1_4_8_sd.dat
    $mycommand -np 32 ./main_2c_19_sd.x 2 4 4 benchw.inp > mio2_19_w_32_2_4_4_sd.dat
    
    echo "Mode: STRONG - Running strong configuration"
    $mycommand -np 32 ./main_1c_19_sd.x 2 4 4 bench1.inp > mio1_19_s_32_2_4_4_sd.dat
    $mycommand -np 32 ./main_1c_19_sd.x 1 4 8 bench1.inp > mio1_19_s_32_1_4_8_sd.dat
    
    $mycommand -np 32 ./main_2c_19_sd.x 2 4 4 bench.inp > mio2_19_s_32_2_4_4_sd.dat
    $mycommand -np 32 ./main_2c_19_sd.x 1 4 8 bench.inp > mio2_19_s_32_1_4_8_sd.dat

fi

if [[ "$mynodes" == "16" ]]; then
    echo "Mode: WEAK - Running weak configuration"
    $mycommand -np 64 ./main_1c_19_sd.x  1 1 64 bench1w.inp > mio1_19_w_64_1_1_64_sd.dat
    $mycommand -np 64 ./main_1c_19_sd.x  1 8 8 bench1w.inp > mio1_19_w_64_1_8_8_sd.dat
    $mycommand -np 64 ./main_1c_19_sd.x  2 4 8 bench1w.inp > mio1_19_w_64_2_4_8_sd.dat
    $mycommand -np 64 ./main_1c_19_sd.x  4 4 4 bench1w.inp > mio1_19_w_64_4_4_4_sd.dat
    
    $mycommand -np 64 ./main_2c_19_sd.x  1 1 64 benchw.inp > mio2_19_w_64_1_1_64_sd.dat
    $mycommand -np 64 ./main_2c_19_sd.x  1 8 8 benchw.inp > mio2_19_w_64_1_8_8_sd.dat
    $mycommand -np 64 ./main_2c_19_sd.x  2 4 8 benchw.inp > mio2_19_w_64_2_4_8_sd.dat
    $mycommand -np 64 ./main_2c_19_sd.x  4 4 4 benchw.inp > mio2_19_w_64_4_4_4_sd.dat
    
    echo "Mode: STRONG - Running strong configuration"
    $mycommand -np 64 ./main_1c_19_sd.x 4 4 4 bench1.inp > mio1_19_s_64_4_4_4_sd.dat
    $mycommand -np 64 ./main_1c_19_sd.x 2 4 8 bench1.inp > mio1_19_s_64_2_4_8_sd.dat
    $mycommand -np 64 ./main_1c_19_sd.x 1 8 8 bench1.inp > mio1_19_s_64_1_8_8_sd.dat
    
    $mycommand -np 64 ./main_2c_19_sd.x 4 4 4 bench.inp > mio2_19_s_64_4_4_4_sd.dat
    $mycommand -np 64 ./main_2c_19_sd.x 2 4 8 bench.inp > mio2_19_s_64_2_4_8_sd.dat
    $mycommand -np 64 ./main_2c_19_sd.x 1 8 8 bench.inp > mio2_19_s_64_1_8_8_sd.dat

fi


if [[ "$mynodes" == "1" ]]; then
    echo "Mode: WEAK - Running weak configuration"
    $mycommand -np 1 ./main_1c_27_sd.x 1 1 1 bench1w.inp > mio1_27_w_1_1_1_1_sd.dat
    $mycommand -np 2 ./main_1c_27_sd.x 1 1 2 bench1w.inp > mio1_27_w_2_1_1_2_sd.dat
    $mycommand -np 4 ./main_1c_27_sd.x 1 2 2 bench1w.inp > mio1_27_w_4_1_2_2_sd.dat
    $mycommand -np 4 ./main_1c_27_sd.x 1 1 4 bench1w.inp > mio1_27_w_4_1_1_4_sd.dat

    $mycommand -np 1 ./main_2c_27_sd.x 1 1 1 benchw.inp > mio2_27_w_1_1_1_1_sd.dat
    $mycommand -np 2 ./main_2c_27_sd.x 1 1 2 benchw.inp > mio2_27_w_2_1_1_2_sd.dat
    $mycommand -np 4 ./main_2c_27_sd.x 1 2 2 benchw.inp > mio2_27_w_4_1_2_2_sd.dat
    $mycommand -np 4 ./main_2c_27_sd.x 1 1 4 benchw.inp > mio2_27_w_4_1_1_4_sd.dat
    
    echo "Mode: STRONG - Running strong configuration"
    
    $mycommand  -np 1 ./main_1c_27_sd.x 1 1 1 bench1.inp > mio1_27_s_1_1_1_1_sd.dat
    $mycommand  -np 2 ./main_1c_27_sd.x 1 1 2 bench1.inp > mio1_27_s_2_1_1_2_sd.dat
    $mycommand  -np 4 ./main_1c_27_sd.x 1 1 4 bench1.inp > mio1_27_s_4_1_1_4_sd.dat
    $mycommand  -np 4 ./main_1c_27_sd.x 1 2 2 bench1.inp > mio1_27_s_4_1_2_2_sd.dat

    $mycommand  -np 1 ./main_2c_27_sd.x 1 1 1 bench.inp > mio2_27_s_1_1_1_1_sd.dat
    $mycommand  -np 2 ./main_2c_27_sd.x 1 1 2 bench.inp > mio2_27_s_2_1_1_2_sd.dat
    $mycommand  -np 4 ./main_2c_27_sd.x 1 1 4 bench.inp > mio2_27_s_4_1_1_4_sd.dat
    $mycommand  -np 4 ./main_2c_27_sd.x 1 2 2 bench.inp > mio2_27_s_4_1_2_2_sd.dat
fi


if [[ "$mynodes" == "2" ]]; then
	echo "Mode: WEAK - Running weak configuration"
        $mycommand -np 8 ./main_1c_27_sd.x 1 1 8 bench1w.inp > mio1_27_w_8_1_1_8_sd.dat
        $mycommand -np 8 ./main_1c_27_sd.x 1 2 4 bench1w.inp > mio1_27_w_8_1_2_4_sd.dat
        $mycommand -np 8 ./main_1c_27_sd.x 2 2 2 bench1w.inp > mio1_27_w_8_2_2_2_sd.dat

	    $mycommand -np 8 ./main_2c_27_sd.x 1 1 8 benchw.inp > mio2_27_w_8_1_1_8_sd.dat
	    $mycommand -np 8 ./main_2c_27_sd.x 1 2 4 benchw.inp > mio2_27_w_8_1_2_4_sd.dat
        $mycommand -np 8 ./main_2c_27_sd.x 2 2 2 benchw.inp > mio2_27_w_8_2_2_2_sd.dat
	    
        echo "Mode: STRONG - Running strong configuration"
        $mycommand -np 8 ./main_1c_27_sd.x  2 2 2 bench1.inp > mio1_27_s_8_2_2_2_sd.dat
        $mycommand -np 8 ./main_1c_27_sd.x  1 2 4 bench1.inp > mio1_27_s_8_1_2_4_sd.dat

        $mycommand -np 8 ./main_2c_27_sd.x  2 2 2 bench.inp > mio2_27_s_8_2_2_2_sd.dat
        $mycommand -np 8 ./main_2c_27_sd.x  1 2 4 bench.inp > mio2_27_s_8_1_2_4_sd.dat

fi

if [[ "$mynodes" == "4" ]]; then
        echo "Mode: WEAK - Running weak configuration"

        $mycommand -np 16 ./main_1c_27_sd.x 1 1 16 bench1w.inp > mio1_27_w_16_1_1_16_sd.dat
        $mycommand -np 16 ./main_1c_27_sd.x 1 4 4 bench1w.inp > mio1_27_w_16_1_4_4_sd.dat
        $mycommand -np 16 ./main_1c_27_sd.x 2 2 4 bench1w.inp > mio1_27_w_16_2_2_4_sd.dat

	    $mycommand -np 16 ./main_2c_27_sd.x 1 1 16 benchw.inp > mio2_27_w_16_1_1_16_sd.dat
	    $mycommand -np 16 ./main_2c_27_sd.x 1 4 4 benchw.inp > mio2_27_w_16_1_4_4_sd.dat
        $mycommand -np 16 ./main_2c_27_sd.x 2 2 4 benchw.inp > mio2_27_w_16_2_2_4_sd.dat
	    
        echo "Mode: STRONG - Running strong configuration"
       $mycommand -np 16 ./main_1c_27_sd.x 2 2 4 bench1.inp > mio1_27_s_16_2_2_4_sd.dat
       $mycommand -np 16 ./main_1c_27_sd.x 1 4 4 bench1.inp > mio1_27_s_16_1_4_4_sd.dat

       $mycommand -np 16 ./main_2c_27_sd.x 2 2 4 bench.inp > mio2_27_s_16_2_2_4_sd.dat
       $mycommand -np 16 ./main_2c_27_sd.x 1 4 4 bench.inp > mio2_27_s_16_1_4_4_sd.dat

fi

if [[ "$mynodes" == "8" ]]; then

    echo "Mode: WEAK - Running weak configuration"
    
    $mycommand -np 32 ./main_1c_27_sd.x 1 1 32 bench1w.inp > mio1_27_w_32_1_1_32_sd.dat
    $mycommand -np 32 ./main_1c_27_sd.x 1 4 8 bench1w.inp > mio1_27_w_32_1_4_8_sd.dat
    $mycommand -np 32 ./main_1c_27_sd.x 2 4 4 bench1w.inp > mio1_27_w_32_2_4_4_sd.dat
    
    $mycommand -np 32 ./main_2c_27_sd.x 1 1 32 benchw.inp > mio2_27_w_32_1_1_32_sd.dat
    $mycommand -np 32 ./main_2c_27_sd.x 1 4 8 benchw.inp > mio2_27_w_32_1_4_8_sd.dat
    $mycommand -np 32 ./main_2c_27_sd.x 2 4 4 benchw.inp > mio2_27_w_32_2_4_4_sd.dat
    
    echo "Mode: STRONG - Running strong configuration"
    $mycommand -np 32 ./main_1c_27_sd.x 2 4 4 bench1.inp > mio1_27_s_32_2_4_4_sd.dat
    $mycommand -np 32 ./main_1c_27_sd.x 1 4 8 bench1.inp > mio1_27_s_32_1_4_8_sd.dat
    
    $mycommand -np 32 ./main_2c_27_sd.x 2 4 4 bench.inp > mio2_27_s_32_2_4_4_sd.dat
    $mycommand -np 32 ./main_2c_27_sd.x 1 4 8 bench.inp > mio2_27_s_32_1_4_8_sd.dat

fi

if [[ "$mynodes" == "16" ]]; then
    echo "Mode: WEAK - Running weak configuration"
    $mycommand -np 64 ./main_1c_27_sd.x  1 1 64 bench1w.inp > mio1_27_w_64_1_1_64_sd.dat
    $mycommand -np 64 ./main_1c_27_sd.x  1 8 8 bench1w.inp > mio1_27_w_64_1_8_8_sd.dat
    $mycommand -np 64 ./main_1c_27_sd.x  2 4 8 bench1w.inp > mio1_27_w_64_2_4_8_sd.dat
    $mycommand -np 64 ./main_1c_27_sd.x  4 4 4 bench1w.inp > mio1_27_w_64_4_4_4_sd.dat
    
    $mycommand -np 64 ./main_2c_27_sd.x  1 1 64 benchw.inp > mio2_27_w_64_1_1_64_sd.dat
    $mycommand -np 64 ./main_2c_27_sd.x  1 8 8 benchw.inp > mio2_27_w_64_1_8_8_sd.dat
    $mycommand -np 64 ./main_2c_27_sd.x  2 4 8 benchw.inp > mio2_27_w_64_2_4_8_sd.dat
    $mycommand -np 64 ./main_2c_27_sd.x  4 4 4 benchw.inp > mio2_27_w_64_4_4_4_sd.dat
    
    echo "Mode: STRONG - Running strong configuration"
    $mycommand -np 64 ./main_1c_27_sd.x 4 4 4 bench1.inp > mio1_27_s_64_4_4_4_sd.dat
    $mycommand -np 64 ./main_1c_27_sd.x 2 4 8 bench1.inp > mio1_27_s_64_2_4_8_sd.dat
    $mycommand -np 64 ./main_1c_27_sd.x 1 8 8 bench1.inp > mio1_27_s_64_1_8_8_sd.dat
    
    $mycommand -np 64 ./main_2c_27_sd.x 4 4 4 bench.inp > mio2_27_s_64_4_4_4_sd.dat
    $mycommand -np 64 ./main_2c_27_sd.x 2 4 8 bench.inp > mio2_27_s_64_2_4_8_sd.dat
    $mycommand -np 64 ./main_2c_27_sd.x 1 8 8 bench.inp > mio2_27_s_64_1_8_8_sd.dat

fi


if [[ "$mynodes" == "1" ]]; then
    echo "Mode: WEAK - Running weak configuration"
    $mycommand -np 1 ./main_1c_27high_sd.x 1 1 1 bench1w.inp > mio1_27h_w_1_1_1_1_sd.dat
    $mycommand -np 2 ./main_1c_27high_sd.x 1 1 2 bench1w.inp > mio1_27h_w_2_1_1_2_sd.dat
    $mycommand -np 4 ./main_1c_27high_sd.x 1 2 2 bench1w.inp > mio1_27h_w_4_1_2_2_sd.dat
    $mycommand -np 4 ./main_1c_27high_sd.x 1 1 4 bench1w.inp > mio1_27h_w_4_1_1_4_sd.dat

    $mycommand -np 1 ./main_2c_27high_sd.x 1 1 1 benchw.inp > mio2_27h_w_1_1_1_1_sd.dat
    $mycommand -np 2 ./main_2c_27high_sd.x 1 1 2 benchw.inp > mio2_27h_w_2_1_1_2_sd.dat
    $mycommand -np 4 ./main_2c_27high_sd.x 1 2 2 benchw.inp > mio2_27h_w_4_1_2_2_sd.dat
    $mycommand -np 4 ./main_2c_27high_sd.x 1 1 4 benchw.inp > mio2_27h_w_4_1_1_4_sd.dat
    
    echo "Mode: STRONG - Running strong configuration"
    
    $mycommand  -np 1 ./main_1c_27high_sd.x 1 1 1 bench1.inp > mio1_27h_s_1_1_1_1_sd.dat
    $mycommand  -np 2 ./main_1c_27high_sd.x 1 1 2 bench1.inp > mio1_27h_s_2_1_1_2_sd.dat
    $mycommand  -np 4 ./main_1c_27high_sd.x 1 1 4 bench1.inp > mio1_27h_s_4_1_1_4_sd.dat
    $mycommand  -np 4 ./main_1c_27high_sd.x 1 2 2 bench1.inp > mio1_27h_s_4_1_2_2_sd.dat

    $mycommand  -np 1 ./main_2c_27high_sd.x 1 1 1 bench.inp > mio2_27h_s_1_1_1_1_sd.dat
    $mycommand  -np 2 ./main_2c_27high_sd.x 1 1 2 bench.inp > mio2_27h_s_2_1_1_2_sd.dat
    $mycommand  -np 4 ./main_2c_27high_sd.x 1 1 4 bench.inp > mio2_27h_s_4_1_1_4_sd.dat
    $mycommand  -np 4 ./main_2c_27high_sd.x 1 2 2 bench.inp > mio2_27h_s_4_1_2_2_sd.dat
fi


if [[ "$mynodes" == "2" ]]; then
	echo "Mode: WEAK - Running weak configuration"
        $mycommand -np 8 ./main_1c_27high_sd.x 1 1 8 bench1w.inp > mio1_27h_w_8_1_1_8_sd.dat
        $mycommand -np 8 ./main_1c_27high_sd.x 1 2 4 bench1w.inp > mio1_27h_w_8_1_2_4_sd.dat
        $mycommand -np 8 ./main_1c_27high_sd.x 2 2 2 bench1w.inp > mio1_27h_w_8_2_2_2_sd.dat

	    $mycommand -np 8 ./main_2c_27high_sd.x 1 1 8 benchw.inp > mio2_27h_w_8_1_1_8_sd.dat
	    $mycommand -np 8 ./main_2c_27high_sd.x 1 2 4 benchw.inp > mio2_27h_w_8_1_2_4_sd.dat
        $mycommand -np 8 ./main_2c_27high_sd.x 2 2 2 benchw.inp > mio2_27h_w_8_2_2_2_sd.dat
	    
        echo "Mode: STRONG - Running strong configuration"
        $mycommand -np 8 ./main_1c_27high_sd.x  2 2 2 bench1.inp > mio1_27h_s_8_2_2_2_sd.dat
        $mycommand -np 8 ./main_1c_27high_sd.x  1 2 4 bench1.inp > mio1_27h_s_8_1_2_4_sd.dat

        $mycommand -np 8 ./main_2c_27high_sd.x  2 2 2 bench.inp > mio2_27h_s_8_2_2_2_sd.dat
        $mycommand -np 8 ./main_2c_27high_sd.x  1 2 4 bench.inp > mio2_27h_s_8_1_2_4_sd.dat

fi

if [[ "$mynodes" == "4" ]]; then
        echo "Mode: WEAK - Running weak configuration"

        $mycommand -np 16 ./main_1c_27high_sd.x 1 1 16 bench1w.inp > mio1_27h_w_16_1_1_16_sd.dat
        $mycommand -np 16 ./main_1c_27high_sd.x 1 4 4 bench1w.inp > mio1_27h_w_16_1_4_4_sd.dat
        $mycommand -np 16 ./main_1c_27high_sd.x 2 2 4 bench1w.inp > mio1_27h_w_16_2_2_4_sd.dat

	    $mycommand -np 16 ./main_2c_27high_sd.x 1 1 16 benchw.inp > mio2_27h_w_16_1_1_16_sd.dat
	    $mycommand -np 16 ./main_2c_27high_sd.x 1 4 4 benchw.inp > mio2_27h_w_16_1_4_4_sd.dat
        $mycommand -np 16 ./main_2c_27high_sd.x 2 2 4 benchw.inp > mio2_27h_w_16_2_2_4_sd.dat
	    
        echo "Mode: STRONG - Running strong configuration"
       $mycommand -np 16 ./main_1c_27high_sd.x 2 2 4 bench1.inp > mio1_27h_s_16_2_2_4_sd.dat
       $mycommand -np 16 ./main_1c_27high_sd.x 1 4 4 bench1.inp > mio1_27h_s_16_1_4_4_sd.dat

       $mycommand -np 16 ./main_2c_27high_sd.x 2 2 4 bench.inp > mio2_27h_s_16_2_2_4_sd.dat
       $mycommand -np 16 ./main_2c_27high_sd.x 1 4 4 bench.inp > mio2_27h_s_16_1_4_4_sd.dat

fi

if [[ "$mynodes" == "8" ]]; then

    echo "Mode: WEAK - Running weak configuration"
    
    $mycommand -np 32 ./main_1c_27high_sd.x 1 1 32 bench1w.inp > mio1_27h_w_32_1_1_32_sd.dat
    $mycommand -np 32 ./main_1c_27high_sd.x 1 4 8 bench1w.inp > mio1_27h_w_32_1_4_8_sd.dat
    $mycommand -np 32 ./main_1c_27high_sd.x 2 4 4 bench1w.inp > mio1_27h_w_32_2_4_4_sd.dat
    
    $mycommand -np 32 ./main_2c_27high_sd.x 1 1 32 benchw.inp > mio2_27h_w_32_1_1_32_sd.dat
    $mycommand -np 32 ./main_2c_27high_sd.x 1 4 8 benchw.inp > mio2_27h_w_32_1_4_8_sd.dat
    $mycommand -np 32 ./main_2c_27high_sd.x 2 4 4 benchw.inp > mio2_27h_w_32_2_4_4_sd.dat
    
    echo "Mode: STRONG - Running strong configuration"
    $mycommand -np 32 ./main_1c_27high_sd.x 2 4 4 bench1.inp > mio1_27h_s_32_2_4_4_sd.dat
    $mycommand -np 32 ./main_1c_27high_sd.x 1 4 8 bench1.inp > mio1_27h_s_32_1_4_8_sd.dat
    
    $mycommand -np 32 ./main_2c_27high_sd.x 2 4 4 bench.inp > mio2_27h_s_32_2_4_4_sd.dat
    $mycommand -np 32 ./main_2c_27high_sd.x 1 4 8 bench.inp > mio2_27h_s_32_1_4_8_sd.dat

fi

if [[ "$mynodes" == "16" ]]; then
    echo "Mode: WEAK - Running weak configuration"
    $mycommand -np 64 ./main_1c_27high_sd.x  1 1 64 bench1w.inp > mio1_27h_w_64_1_1_64_sd.dat
    $mycommand -np 64 ./main_1c_27high_sd.x  1 8 8 bench1w.inp > mio1_27h_w_64_1_8_8_sd.dat
    $mycommand -np 64 ./main_1c_27high_sd.x  2 4 8 bench1w.inp > mio1_27h_w_64_2_4_8_sd.dat
    $mycommand -np 64 ./main_1c_27high_sd.x  4 4 4 bench1w.inp > mio1_27h_w_64_4_4_4_sd.dat
    
    $mycommand -np 64 ./main_2c_27high_sd.x  1 1 64 benchw.inp > mio2_27h_w_64_1_1_64_sd.dat
    $mycommand -np 64 ./main_2c_27high_sd.x  1 8 8 benchw.inp > mio2_27h_w_64_1_8_8_sd.dat
    $mycommand -np 64 ./main_2c_27high_sd.x  2 4 8 benchw.inp > mio2_27h_w_64_2_4_8_sd.dat
    $mycommand -np 64 ./main_2c_27high_sd.x  4 4 4 benchw.inp > mio2_27h_w_64_4_4_4_sd.dat
    
    echo "Mode: STRONG - Running strong configuration"
    $mycommand -np 64 ./main_1c_27high_sd.x 4 4 4 bench1.inp > mio1_27h_s_64_4_4_4_sd.dat
    $mycommand -np 64 ./main_1c_27high_sd.x 2 4 8 bench1.inp > mio1_27h_s_64_2_4_8_sd.dat
    $mycommand -np 64 ./main_1c_27high_sd.x 1 8 8 bench1.inp > mio1_27h_s_64_1_8_8_sd.dat
    
    $mycommand -np 64 ./main_2c_27high_sd.x 4 4 4 bench.inp > mio2_27h_s_64_4_4_4_sd.dat
    $mycommand -np 64 ./main_2c_27high_sd.x 2 4 8 bench.inp > mio2_27h_s_64_2_4_8_sd.dat
    $mycommand -np 64 ./main_2c_27high_sd.x 1 8 8 bench.inp > mio2_27h_s_64_1_8_8_sd.dat

fi


if [[ "$mynodes" == "32" ]]; then
        echo "Mode: WEAK - Running weak configuration"

        $mycommand -np 128 ./main_1c_15.x 1 1 128 bench1w.inp > mio1_15_w_128_1_1_128.dat
        $mycommand -np 128 ./main_1c_15.x 1 8 16 bench1w.inp > mio1_15_w_128_1_8_16.dat
        $mycommand -np 128 ./main_1c_15.x 4 4 8 bench1w.inp > mio1_15_w_128_4_4_8.dat

	    $mycommand -np 128 ./main_2c_15.x 1 1 128 benchw.inp > mio2_15_w_128_1_1_128.dat
	    $mycommand -np 128 ./main_2c_15.x 1 8 16 benchw.inp > mio2_15_w_128_1_8_16.dat
        $mycommand -np 128 ./main_2c_15.x 4 4 8 benchw.inp > mio2_15_w_128_4_4_8.dat
	    
        echo "Mode: STRONG - Running strong configuration"
       $mycommand -np 128 ./main_1c_15.x 4 4 8 bench1.inp > mio1_15_s_128_4_4_8.dat
       $mycommand -np 128 ./main_1c_15.x 1 8 16 bench1.inp > mio1_15_s_128_1_8_16.dat

       $mycommand -np 128 ./main_2c_15.x 4 4 8 bench.inp > mio2_15_s_128_4_4_8.dat
       $mycommand -np 128 ./main_2c_15.x 1 8 16 bench.inp > mio2_15_s_128_1_8_16.dat

fi

if [[ "$mynodes" == "64" ]]; then

    echo "Mode: WEAK - Running weak configuration"
    
    $mycommand -np 256 ./main_1c_15.x 1 1 256 bench1w.inp > mio1_15_w_256_1_1_256.dat
    $mycommand -np 256 ./main_1c_15.x 1 16 16 bench1w.inp > mio1_15_w_256_1_16_16.dat
    $mycommand -np 256 ./main_1c_15.x 4 8 8 bench1w.inp > mio1_15_w_256_4_8_8.dat
    
    $mycommand -np 256 ./main_2c_15.x 1 1 256 benchw.inp > mio2_15_w_256_1_1_256.dat
    $mycommand -np 256 ./main_2c_15.x 1 16 16 benchw.inp > mio2_15_w_256_1_16_16.dat
    $mycommand -np 256 ./main_2c_15.x 4 8 8 benchw.inp > mio2_15_w_256_4_8_8.dat
    
    echo "Mode: STRONG - Running strong configuration"
    $mycommand -np 256 ./main_1c_15.x 4 8 8 bench1.inp > mio1_15_s_256_4_8_8.dat
    $mycommand -np 256 ./main_1c_15.x 1 16 16 bench1.inp > mio1_15_s_256_1_16_16.dat
    
    $mycommand -np 256 ./main_2c_15.x 4 8 8 bench.inp > mio2_15_s_256_4_8_8.dat
    $mycommand -np 256 ./main_2c_15.x 1 16 16 bench.inp > mio2_15_s_256_1_16_16.dat

fi

if [[ "$mynodes" == "128" ]]; then
    echo "Mode: WEAK - Running weak configuration"
    $mycommand -np 512 ./main_1c_15.x  1 1 512 bench1w.inp > mio1_15_w_512_1_1_512.dat
    $mycommand -np 512 ./main_1c_15.x  1 16 32 bench1w.inp > mio1_15_w_512_1_16_32.dat
    $mycommand -np 512 ./main_1c_15.x  4 8 16 bench1w.inp > mio1_15_w_512_4_8_16.dat
    $mycommand -np 512 ./main_1c_15.x  8 8 8 bench1w.inp > mio1_15_w_512_8_8_8.dat
    
    $mycommand -np 512 ./main_2c_15.x  1 1 512 benchw.inp > mio2_15_w_512_1_1_512.dat
    $mycommand -np 512 ./main_2c_15.x  1 16 32 benchw.inp > mio2_15_w_512_1_16_32.dat
    $mycommand -np 512 ./main_2c_15.x  4 8 16 benchw.inp > mio2_15_w_512_4_8_16.dat
    $mycommand -np 512 ./main_2c_15.x  8 8 8 benchw.inp > mio2_15_w_512_8_8_8.dat
    
    echo "Mode: STRONG - Running strong configuration"
    $mycommand -np 512 ./main_1c_15.x 8 8 8 bench1.inp > mio1_15_s_512_8_8_8.dat
    $mycommand -np 512 ./main_1c_15.x 4 8 16 bench1.inp > mio1_15_s_512_4_8_16.dat
    $mycommand -np 512 ./main_1c_15.x 1 16 32 bench1.inp > mio1_15_s_512_1_16_32.dat
    
    $mycommand -np 512 ./main_2c_15.x 8 8 8 bench.inp > mio2_15_s_512_8_8_8.dat
    $mycommand -np 512 ./main_2c_15.x 4 8 16 bench.inp > mio2_15_s_512_4_8_16.dat
    $mycommand -np 512 ./main_2c_15.x 1 16 32 bench.inp > mio2_15_s_512_1_16_32.dat

fi

if [[ "$mynodes" == "32" ]]; then
        echo "Mode: WEAK - Running weak configuration"
       
        $mycommand -np 128 ./main_1c_19.x 1 1 128 bench1w.inp > mio1_19_w_128_1_1_128.dat
        $mycommand -np 128 ./main_1c_19.x 1 8 16 bench1w.inp > mio1_19_w_128_1_8_16.dat
        $mycommand -np 128 ./main_1c_19.x 4 4 8 bench1w.inp > mio1_19_w_128_4_4_8.dat

	    $mycommand -np 128 ./main_2c_19.x 1 1 128 benchw.inp > mio2_19_w_128_1_1_128.dat
	    $mycommand -np 128 ./main_2c_19.x 1 8 16 benchw.inp > mio2_19_w_128_1_8_16.dat
        $mycommand -np 128 ./main_2c_19.x 4 4 8 benchw.inp > mio2_19_w_128_4_4_8.dat
	    
        echo "Mode: STRONG - Running strong configuration"
       $mycommand -np 128 ./main_1c_19.x 4 4 8 bench1.inp > mio1_19_s_128_4_4_8.dat
       $mycommand -np 128 ./main_1c_19.x 1 8 16 bench1.inp > mio1_19_s_128_1_8_16.dat

       $mycommand -np 128 ./main_2c_19.x 4 4 8 bench.inp > mio2_19_s_128_4_4_8.dat
       $mycommand -np 128 ./main_2c_19.x 1 8 16 bench.inp > mio2_19_s_128_1_8_16.dat

fi

if [[ "$mynodes" == "64" ]]; then

    echo "Mode: WEAK - Running weak configuration"
    
    $mycommand -np 256 ./main_1c_19.x 1 1 256 bench1w.inp > mio1_19_w_256_1_1_256.dat
    $mycommand -np 256 ./main_1c_19.x 1 16 16 bench1w.inp > mio1_19_w_256_1_16_16.dat
    $mycommand -np 256 ./main_1c_19.x 4 8 8 bench1w.inp > mio1_19_w_256_4_8_8.dat
    
    $mycommand -np 256 ./main_2c_19.x 1 1 256 benchw.inp > mio2_19_w_256_1_1_256.dat
    $mycommand -np 256 ./main_2c_19.x 1 16 16 benchw.inp > mio2_19_w_256_1_16_16.dat
    $mycommand -np 256 ./main_2c_19.x 4 8 8 benchw.inp > mio2_19_w_256_4_8_8.dat
    
    echo "Mode: STRONG - Running strong configuration"
    $mycommand -np 256 ./main_1c_19.x 4 8 8 bench1.inp > mio1_19_s_256_4_8_8.dat
    $mycommand -np 256 ./main_1c_19.x 1 16 16 bench1.inp > mio1_19_s_256_1_16_16.dat
    
    $mycommand -np 256 ./main_2c_19.x 4 8 8 bench.inp > mio2_19_s_256_4_8_8.dat
    $mycommand -np 256 ./main_2c_19.x 1 16 16 bench.inp > mio2_19_s_256_1_16_16.dat

fi

if [[ "$mynodes" == "128" ]]; then
    echo "Mode: WEAK - Running weak configuration"
    
    $mycommand -np 512 ./main_1c_19.x  1 1 512 bench1w.inp > mio1_19_w_512_1_1_512.dat
    $mycommand -np 512 ./main_1c_19.x  1 16 32 bench1w.inp > mio1_19_w_512_1_16_32.dat
    $mycommand -np 512 ./main_1c_19.x  4 8 16 bench1w.inp > mio1_19_w_512_4_8_16.dat
    $mycommand -np 512 ./main_1c_19.x  8 8 8 bench1w.inp > mio1_19_w_512_8_8_8.dat
    
    $mycommand -np 512 ./main_2c_19.x  1 1 512 benchw.inp > mio2_19_w_512_1_1_512.dat
    $mycommand -np 512 ./main_2c_19.x  1 16 32 benchw.inp > mio2_19_w_512_1_16_32.dat
    $mycommand -np 512 ./main_2c_19.x  4 8 16 benchw.inp > mio2_19_w_512_4_8_16.dat
    $mycommand -np 512 ./main_2c_19.x  8 8 8 benchw.inp > mio2_19_w_512_8_8_8.dat
    
    echo "Mode: STRONG - Running strong configuration"
    $mycommand -np 512 ./main_1c_19.x 8 8 8 bench1.inp > mio1_19_s_512_8_8_8.dat
    $mycommand -np 512 ./main_1c_19.x 4 8 16 bench1.inp > mio1_19_s_512_4_8_16.dat
    $mycommand -np 512 ./main_1c_19.x 1 16 32 bench1.inp > mio1_19_s_512_1_16_32.dat
    
    $mycommand -np 512 ./main_2c_19.x 8 8 8 bench.inp > mio2_19_s_512_8_8_8.dat
    $mycommand -np 512 ./main_2c_19.x 4 8 16 bench.inp > mio2_19_s_512_4_8_16.dat
    $mycommand -np 512 ./main_2c_19.x 1 16 32 bench.inp > mio2_19_s_512_1_16_32.dat


fi


if [[ "$mynodes" == "32" ]]; then
        echo "Mode: WEAK - Running weak configuration"

        $mycommand -np 128 ./main_1c_27.x 1 1 128 bench1w.inp > mio1_27_w_128_1_1_128.dat
        $mycommand -np 128 ./main_1c_27.x 1 8 16 bench1w.inp > mio1_27_w_128_1_8_16.dat
        $mycommand -np 128 ./main_1c_27.x 4 4 8 bench1w.inp > mio1_27_w_128_4_4_8.dat

	    $mycommand -np 128 ./main_2c_27.x 1 1 128 benchw.inp > mio2_27_w_128_1_1_128.dat
	    $mycommand -np 128 ./main_2c_27.x 1 8 16 benchw.inp > mio2_27_w_128_1_8_16.dat
        $mycommand -np 128 ./main_2c_27.x 4 4 8 benchw.inp > mio2_27_w_128_4_4_8.dat
	    
        echo "Mode: STRONG - Running strong configuration"
       $mycommand -np 128 ./main_1c_27.x 4 4 8 bench1.inp > mio1_27_s_128_4_4_8.dat
       $mycommand -np 128 ./main_1c_27.x 1 8 16 bench1.inp > mio1_27_s_128_1_8_16.dat

       $mycommand -np 128 ./main_2c_27.x 4 4 8 bench.inp > mio2_27_s_128_4_4_8.dat
       $mycommand -np 128 ./main_2c_27.x 1 8 16 bench.inp > mio2_27_s_128_1_8_16.dat

fi

if [[ "$mynodes" == "64" ]]; then

    echo "Mode: WEAK - Running weak configuration"
    
    $mycommand -np 256 ./main_1c_27.x 1 1 256 bench1w.inp > mio1_27_w_256_1_1_256.dat
    $mycommand -np 256 ./main_1c_27.x 1 16 16 bench1w.inp > mio1_27_w_256_1_16_16.dat
    $mycommand -np 256 ./main_1c_27.x 4 8 8 bench1w.inp > mio1_27_w_256_4_8_8.dat
    
    $mycommand -np 256 ./main_2c_27.x 1 1 256 benchw.inp > mio2_27_w_256_1_1_256.dat
    $mycommand -np 256 ./main_2c_27.x 1 16 16 benchw.inp > mio2_27_w_256_1_16_16.dat
    $mycommand -np 256 ./main_2c_27.x 4 8 8 benchw.inp > mio2_27_w_256_4_8_8.dat
    
    echo "Mode: STRONG - Running strong configuration"
    $mycommand -np 256 ./main_1c_27.x 4 8 8 bench1.inp > mio1_27_s_256_4_8_8.dat
    $mycommand -np 256 ./main_1c_27.x 1 16 16 bench1.inp > mio1_27_s_256_1_16_16.dat
    
    $mycommand -np 256 ./main_2c_27.x 4 8 8 bench.inp > mio2_27_s_256_4_8_8.dat
    $mycommand -np 256 ./main_2c_27.x 1 16 16 bench.inp > mio2_27_s_256_1_16_16.dat

fi

if [[ "$mynodes" == "128" ]]; then
    echo "Mode: WEAK - Running weak configuration"

    $mycommand -np 512 ./main_1c_27.x  1 1 512 bench1w.inp > mio1_27_w_512_1_1_512.dat
    $mycommand -np 512 ./main_1c_27.x  1 16 32 bench1w.inp > mio1_27_w_512_1_16_32.dat
    $mycommand -np 512 ./main_1c_27.x  4 8 16 bench1w.inp > mio1_27_w_512_4_8_16.dat
    $mycommand -np 512 ./main_1c_27.x  8 8 8 bench1w.inp > mio1_27_w_512_8_8_8.dat
    
    $mycommand -np 512 ./main_2c_27.x  1 1 512 benchw.inp > mio2_27_w_512_1_1_512.dat
    $mycommand -np 512 ./main_2c_27.x  1 16 32 benchw.inp > mio2_27_w_512_1_16_32.dat
    $mycommand -np 512 ./main_2c_27.x  4 8 16 benchw.inp > mio2_27_w_512_4_8_16.dat
    $mycommand -np 512 ./main_2c_27.x  8 8 8 benchw.inp > mio2_27_w_512_8_8_8.dat
    
    echo "Mode: STRONG - Running strong configuration"
    $mycommand -np 512 ./main_1c_27.x 8 8 8 bench1.inp > mio1_27_s_512_8_8_8.dat
    $mycommand -np 512 ./main_1c_27.x 4 8 16 bench1.inp > mio1_27_s_512_4_8_16.dat
    $mycommand -np 512 ./main_1c_27.x 1 16 32 bench1.inp > mio1_27_s_512_1_16_32.dat
    
    $mycommand -np 512 ./main_2c_27.x 8 8 8 bench.inp > mio2_27_s_512_8_8_8.dat
    $mycommand -np 512 ./main_2c_27.x 4 8 16 bench.inp > mio2_27_s_512_4_8_16.dat
    $mycommand -np 512 ./main_2c_27.x 1 16 32 bench.inp > mio2_27_s_512_1_16_32.dat

fi


if [[ "$mynodes" == "32" ]]; then
        echo "Mode: WEAK - Running weak configuration"

        $mycommand -np 128 ./main_1c_27high.x 1 1 128 bench1w.inp > mio1_27h_w_128_1_1_128.dat
        $mycommand -np 128 ./main_1c_27high.x 1 8 16 bench1w.inp > mio1_27h_w_128_1_8_16.dat
        $mycommand -np 128 ./main_1c_27high.x 4 4 8 bench1w.inp > mio1_27h_w_128_4_4_8.dat

	    $mycommand -np 128 ./main_2c_27high.x 1 1 128 benchw.inp > mio2_27h_w_128_1_1_128.dat
	    $mycommand -np 128 ./main_2c_27high.x 1 8 16 benchw.inp > mio2_27h_w_128_1_8_16.dat
        $mycommand -np 128 ./main_2c_27high.x 4 4 8 benchw.inp > mio2_27h_w_128_4_4_8.dat
	    
        echo "Mode: STRONG - Running strong configuration"
       $mycommand -np 128 ./main_1c_27high.x 4 4 8 bench1.inp > mio1_27h_s_128_4_4_8.dat
       $mycommand -np 128 ./main_1c_27high.x 1 8 16 bench1.inp > mio1_27h_s_128_1_8_16.dat

       $mycommand -np 128 ./main_2c_27high.x 4 4 8 bench.inp > mio2_27h_s_128_4_4_8.dat
       $mycommand -np 128 ./main_2c_27high.x 1 8 16 bench.inp > mio2_27h_s_128_1_8_16.dat

fi

if [[ "$mynodes" == "64" ]]; then

    echo "Mode: WEAK - Running weak configuration"
    
    $mycommand -np 256 ./main_1c_27high.x 1 1 256 bench1w.inp > mio1_27h_w_256_1_1_256.dat
    $mycommand -np 256 ./main_1c_27high.x 1 16 16 bench1w.inp > mio1_27h_w_256_1_16_16.dat
    $mycommand -np 256 ./main_1c_27high.x 4 8 8 bench1w.inp > mio1_27h_w_256_4_8_8.dat
    
    $mycommand -np 256 ./main_2c_27high.x 1 1 256 benchw.inp > mio2_27h_w_256_1_1_256.dat
    $mycommand -np 256 ./main_2c_27high.x 1 16 16 benchw.inp > mio2_27h_w_256_1_16_16.dat
    $mycommand -np 256 ./main_2c_27high.x 4 8 8 benchw.inp > mio2_27h_w_256_4_8_8.dat
    
    echo "Mode: STRONG - Running strong configuration"
    $mycommand -np 256 ./main_1c_27high.x 4 8 8 bench1.inp > mio1_27h_s_256_4_8_8.dat
    $mycommand -np 256 ./main_1c_27high.x 1 16 16 bench1.inp > mio1_27h_s_256_1_16_16.dat
    
    $mycommand -np 256 ./main_2c_27high.x 4 8 8 bench.inp > mio2_27h_s_256_4_8_8.dat
    $mycommand -np 256 ./main_2c_27high.x 1 16 16 bench.inp > mio2_27h_s_256_1_16_16.dat

fi

if [[ "$mynodes" == "128" ]]; then
    echo "Mode: WEAK - Running weak configuration"

    $mycommand -np 512 ./main_1c_27high.x  1 1 512 bench1w.inp > mio1_27h_w_512_1_1_512.dat
    $mycommand -np 512 ./main_1c_27high.x  1 16 32 bench1w.inp > mio1_27h_w_512_1_16_32.dat
    $mycommand -np 512 ./main_1c_27high.x  4 8 16 bench1w.inp > mio1_27h_w_512_4_8_16.dat
    $mycommand -np 512 ./main_1c_27high.x  8 8 8 bench1w.inp > mio1_27h_w_512_8_8_8.dat
    
    $mycommand -np 512 ./main_2c_27high.x  1 1 512 benchw.inp > mio2_27h_w_512_1_1_512.dat
    $mycommand -np 512 ./main_2c_27high.x  1 16 32 benchw.inp > mio2_27h_w_512_1_16_32.dat
    $mycommand -np 512 ./main_2c_27high.x  4 8 16 benchw.inp > mio2_27h_w_512_4_8_16.dat
    $mycommand -np 512 ./main_2c_27high.x  8 8 8 benchw.inp > mio2_27h_w_512_8_8_8.dat
    
    echo "Mode: STRONG - Running strong configuration"
    $mycommand -np 512 ./main_1c_27high.x 8 8 8 bench1.inp > mio1_27h_s_512_8_8_8.dat
    $mycommand -np 512 ./main_1c_27high.x 4 8 16 bench1.inp > mio1_27h_s_512_4_8_16.dat
    $mycommand -np 512 ./main_1c_27high.x 1 16 32 bench1.inp > mio1_27h_s_512_1_16_32.dat
    
    $mycommand -np 512 ./main_2c_27high.x 8 8 8 bench.inp > mio2_27h_s_512_8_8_8.dat
    $mycommand -np 512 ./main_2c_27high.x 4 8 16 bench.inp > mio2_27h_s_512_4_8_16.dat
    $mycommand -np 512 ./main_2c_27high.x 1 16 32 bench.inp > mio2_27h_s_512_1_16_32.dat

fi

##########################################################################################################################

if [[ "$mynodes" == "32" ]]; then
        echo "Mode: WEAK - Running weak configuration"

        $mycommand -np 128 ./main_1c_15_d.x 1 1 128 bench1w.inp > mio1_15_w_128_1_1_128_d.dat
        $mycommand -np 128 ./main_1c_15_d.x 1 8 16 bench1w.inp > mio1_15_w_128_1_8_16_d.dat
        $mycommand -np 128 ./main_1c_15_d.x 4 4 8 bench1w.inp > mio1_15_w_128_4_4_8_d.dat

	    $mycommand -np 128 ./main_2c_15_d.x 1 1 128 benchw.inp > mio2_15_w_128_1_1_128_d.dat
	    $mycommand -np 128 ./main_2c_15_d.x 1 8 16 benchw.inp > mio2_15_w_128_1_8_16_d.dat
        $mycommand -np 128 ./main_2c_15_d.x 4 4 8 benchw.inp > mio2_15_w_128_4_4_8_d.dat
	    
        echo "Mode: STRONG - Running strong configuration"
       $mycommand -np 128 ./main_1c_15_d.x 4 4 8 bench1.inp > mio1_15_s_128_4_4_8_d.dat
       $mycommand -np 128 ./main_1c_15_d.x 1 8 16 bench1.inp > mio1_15_s_128_1_8_16_d.dat

       $mycommand -np 128 ./main_2c_15_d.x 4 4 8 bench.inp > mio2_15_s_128_4_4_8_d.dat
       $mycommand -np 128 ./main_2c_15_d.x 1 8 16 bench.inp > mio2_15_s_128_1_8_16_d.dat

fi

if [[ "$mynodes" == "64" ]]; then

    echo "Mode: WEAK - Running weak configuration"
    
    $mycommand -np 256 ./main_1c_15_d.x 1 1 256 bench1w.inp > mio1_15_w_256_1_1_256_d.dat
    $mycommand -np 256 ./main_1c_15_d.x 1 16 16 bench1w.inp > mio1_15_w_256_1_16_16_d.dat
    $mycommand -np 256 ./main_1c_15_d.x 4 8 8 bench1w.inp > mio1_15_w_256_4_8_8_d.dat
    
    $mycommand -np 256 ./main_2c_15_d.x 1 1 256 benchw.inp > mio2_15_w_256_1_1_256_d.dat
    $mycommand -np 256 ./main_2c_15_d.x 1 16 16 benchw.inp > mio2_15_w_256_1_16_16_d.dat
    $mycommand -np 256 ./main_2c_15_d.x 4 8 8 benchw.inp > mio2_15_w_256_4_8_8_d.dat
    
    echo "Mode: STRONG - Running strong configuration"
    $mycommand -np 256 ./main_1c_15_d.x 4 8 8 bench1.inp > mio1_15_s_256_4_8_8_d.dat
    $mycommand -np 256 ./main_1c_15_d.x 1 16 16 bench1.inp > mio1_15_s_256_1_16_16_d.dat
    
    $mycommand -np 256 ./main_2c_15_d.x 4 8 8 bench.inp > mio2_15_s_256_4_8_8_d.dat
    $mycommand -np 256 ./main_2c_15_d.x 1 16 16 bench.inp > mio2_15_s_256_1_16_16_d.dat

fi

if [[ "$mynodes" == "128" ]]; then
    echo "Mode: WEAK - Running weak configuration"
    $mycommand -np 512 ./main_1c_15_d.x  1 1 512 bench1w.inp > mio1_15_w_512_1_1_512_d.dat
    $mycommand -np 512 ./main_1c_15_d.x  1 16 32 bench1w.inp > mio1_15_w_512_1_16_32_d.dat
    $mycommand -np 512 ./main_1c_15_d.x  4 8 16 bench1w.inp > mio1_15_w_512_4_8_16_d.dat
    $mycommand -np 512 ./main_1c_15_d.x  8 8 8 bench1w.inp > mio1_15_w_512_8_8_8_d.dat
    
    $mycommand -np 512 ./main_2c_15_d.x  1 1 512 benchw.inp > mio2_15_w_512_1_1_512_d.dat
    $mycommand -np 512 ./main_2c_15_d.x  1 16 32 benchw.inp > mio2_15_w_512_1_16_32_d.dat
    $mycommand -np 512 ./main_2c_15_d.x  4 8 16 benchw.inp > mio2_15_w_512_4_8_16_d.dat
    $mycommand -np 512 ./main_2c_15_d.x  8 8 8 benchw.inp > mio2_15_w_512_8_8_8_d.dat
    
    echo "Mode: STRONG - Running strong configuration"
    $mycommand -np 512 ./main_1c_15_d.x 8 8 8 bench1.inp > mio1_15_s_512_8_8_8_d.dat
    $mycommand -np 512 ./main_1c_15_d.x 4 8 16 bench1.inp > mio1_15_s_512_4_8_16_d.dat
    $mycommand -np 512 ./main_1c_15_d.x 1 16 32 bench1.inp > mio1_15_s_512_1_16_32_d.dat
    
    $mycommand -np 512 ./main_2c_15_d.x 8 8 8 bench.inp > mio2_15_s_512_8_8_8_d.dat
    $mycommand -np 512 ./main_2c_15_d.x 4 8 16 bench.inp > mio2_15_s_512_4_8_16_d.dat
    $mycommand -np 512 ./main_2c_15_d.x 1 16 32 bench.inp > mio2_15_s_512_1_16_32_d.dat

fi

if [[ "$mynodes" == "32" ]]; then
        echo "Mode: WEAK - Running weak configuration"
       
        $mycommand -np 128 ./main_1c_19_d.x 1 1 128 bench1w.inp > mio1_19_w_128_1_1_128_d.dat
        $mycommand -np 128 ./main_1c_19_d.x 1 8 16 bench1w.inp > mio1_19_w_128_1_8_16_d.dat
        $mycommand -np 128 ./main_1c_19_d.x 4 4 8 bench1w.inp > mio1_19_w_128_4_4_8_d.dat

	    $mycommand -np 128 ./main_2c_19_d.x 1 1 128 benchw.inp > mio2_19_w_128_1_1_128_d.dat
	    $mycommand -np 128 ./main_2c_19_d.x 1 8 16 benchw.inp > mio2_19_w_128_1_8_16_d.dat
        $mycommand -np 128 ./main_2c_19_d.x 4 4 8 benchw.inp > mio2_19_w_128_4_4_8_d.dat
	    
        echo "Mode: STRONG - Running strong configuration"
       $mycommand -np 128 ./main_1c_19_d.x 4 4 8 bench1.inp > mio1_19_s_128_4_4_8_d.dat
       $mycommand -np 128 ./main_1c_19_d.x 1 8 16 bench1.inp > mio1_19_s_128_1_8_16_d.dat

       $mycommand -np 128 ./main_2c_19_d.x 4 4 8 bench.inp > mio2_19_s_128_4_4_8_d.dat
       $mycommand -np 128 ./main_2c_19_d.x 1 8 16 bench.inp > mio2_19_s_128_1_8_16_d.dat

fi

if [[ "$mynodes" == "64" ]]; then

    echo "Mode: WEAK - Running weak configuration"
    
    $mycommand -np 256 ./main_1c_19_d.x 1 1 256 bench1w.inp > mio1_19_w_256_1_1_256_d.dat
    $mycommand -np 256 ./main_1c_19_d.x 1 16 16 bench1w.inp > mio1_19_w_256_1_16_16_d.dat
    $mycommand -np 256 ./main_1c_19_d.x 4 8 8 bench1w.inp > mio1_19_w_256_4_8_8_d.dat
    
    $mycommand -np 256 ./main_2c_19_d.x 1 1 256 benchw.inp > mio2_19_w_256_1_1_256_d.dat
    $mycommand -np 256 ./main_2c_19_d.x 1 16 16 benchw.inp > mio2_19_w_256_1_16_16_d.dat
    $mycommand -np 256 ./main_2c_19_d.x 4 8 8 benchw.inp > mio2_19_w_256_4_8_8_d.dat
    
    echo "Mode: STRONG - Running strong configuration"
    $mycommand -np 256 ./main_1c_19_d.x 4 8 8 bench1.inp > mio1_19_s_256_4_8_8_d.dat
    $mycommand -np 256 ./main_1c_19_d.x 1 16 16 bench1.inp > mio1_19_s_256_1_16_16_d.dat
    
    $mycommand -np 256 ./main_2c_19_d.x 4 8 8 bench.inp > mio2_19_s_256_4_8_8_d.dat
    $mycommand -np 256 ./main_2c_19_d.x 1 16 16 bench.inp > mio2_19_s_256_1_16_16_d.dat

fi

if [[ "$mynodes" == "128" ]]; then
    echo "Mode: WEAK - Running weak configuration"
    
    $mycommand -np 512 ./main_1c_19_d.x  1 1 512 bench1w.inp > mio1_19_w_512_1_1_512_d.dat
    $mycommand -np 512 ./main_1c_19_d.x  1 16 32 bench1w.inp > mio1_19_w_512_1_16_32_d.dat
    $mycommand -np 512 ./main_1c_19_d.x  4 8 16 bench1w.inp > mio1_19_w_512_4_8_16_d.dat
    $mycommand -np 512 ./main_1c_19_d.x  8 8 8 bench1w.inp > mio1_19_w_512_8_8_8_d.dat
    
    $mycommand -np 512 ./main_2c_19_d.x  1 1 512 benchw.inp > mio2_19_w_512_1_1_512_d.dat
    $mycommand -np 512 ./main_2c_19_d.x  1 16 32 benchw.inp > mio2_19_w_512_1_16_32_d.dat
    $mycommand -np 512 ./main_2c_19_d.x  4 8 16 benchw.inp > mio2_19_w_512_4_8_16_d.dat
    $mycommand -np 512 ./main_2c_19_d.x  8 8 8 benchw.inp > mio2_19_w_512_8_8_8_d.dat
    
    echo "Mode: STRONG - Running strong configuration"
    $mycommand -np 512 ./main_1c_19_d.x 8 8 8 bench1.inp > mio1_19_s_512_8_8_8_d.dat
    $mycommand -np 512 ./main_1c_19_d.x 4 8 16 bench1.inp > mio1_19_s_512_4_8_16_d.dat
    $mycommand -np 512 ./main_1c_19_d.x 1 16 32 bench1.inp > mio1_19_s_512_1_16_32_d.dat
    
    $mycommand -np 512 ./main_2c_19_d.x 8 8 8 bench.inp > mio2_19_s_512_8_8_8_d.dat
    $mycommand -np 512 ./main_2c_19_d.x 4 8 16 bench.inp > mio2_19_s_512_4_8_16_d.dat
    $mycommand -np 512 ./main_2c_19_d.x 1 16 32 bench.inp > mio2_19_s_512_1_16_32_d.dat


fi


if [[ "$mynodes" == "32" ]]; then
        echo "Mode: WEAK - Running weak configuration"

        $mycommand -np 128 ./main_1c_27_d.x 1 1 128 bench1w.inp > mio1_27_w_128_1_1_128_d.dat
        $mycommand -np 128 ./main_1c_27_d.x 1 8 16 bench1w.inp > mio1_27_w_128_1_8_16_d.dat
        $mycommand -np 128 ./main_1c_27_d.x 4 4 8 bench1w.inp > mio1_27_w_128_4_4_8_d.dat

	    $mycommand -np 128 ./main_2c_27_d.x 1 1 128 benchw.inp > mio2_27_w_128_1_1_128_d.dat
	    $mycommand -np 128 ./main_2c_27_d.x 1 8 16 benchw.inp > mio2_27_w_128_1_8_16_d.dat
        $mycommand -np 128 ./main_2c_27_d.x 4 4 8 benchw.inp > mio2_27_w_128_4_4_8_d.dat
	    
        echo "Mode: STRONG - Running strong configuration"
       $mycommand -np 128 ./main_1c_27_d.x 4 4 8 bench1.inp > mio1_27_s_128_4_4_8_d.dat
       $mycommand -np 128 ./main_1c_27_d.x 1 8 16 bench1.inp > mio1_27_s_128_1_8_16_d.dat

       $mycommand -np 128 ./main_2c_27_d.x 4 4 8 bench.inp > mio2_27_s_128_4_4_8_d.dat
       $mycommand -np 128 ./main_2c_27_d.x 1 8 16 bench.inp > mio2_27_s_128_1_8_16_d.dat

fi

if [[ "$mynodes" == "64" ]]; then

    echo "Mode: WEAK - Running weak configuration"
    
    $mycommand -np 256 ./main_1c_27_d.x 1 1 256 bench1w.inp > mio1_27_w_256_1_1_256_d.dat
    $mycommand -np 256 ./main_1c_27_d.x 1 16 16 bench1w.inp > mio1_27_w_256_1_16_16_d.dat
    $mycommand -np 256 ./main_1c_27_d.x 4 8 8 bench1w.inp > mio1_27_w_256_4_8_8_d.dat
    
    $mycommand -np 256 ./main_2c_27_d.x 1 1 256 benchw.inp > mio2_27_w_256_1_1_256_d.dat
    $mycommand -np 256 ./main_2c_27_d.x 1 16 16 benchw.inp > mio2_27_w_256_1_16_16_d.dat
    $mycommand -np 256 ./main_2c_27_d.x 4 8 8 benchw.inp > mio2_27_w_256_4_8_8_d.dat
    
    echo "Mode: STRONG - Running strong configuration"
    $mycommand -np 256 ./main_1c_27_d.x 4 8 8 bench1.inp > mio1_27_s_256_4_8_8_d.dat
    $mycommand -np 256 ./main_1c_27_d.x 1 16 16 bench1.inp > mio1_27_s_256_1_16_16_d.dat
    
    $mycommand -np 256 ./main_2c_27_d.x 4 8 8 bench.inp > mio2_27_s_256_4_8_8_d.dat
    $mycommand -np 256 ./main_2c_27_d.x 1 16 16 bench.inp > mio2_27_s_256_1_16_16_d.dat

fi

if [[ "$mynodes" == "128" ]]; then
    echo "Mode: WEAK - Running weak configuration"

    $mycommand -np 512 ./main_1c_27_d.x  1 1 512 bench1w.inp > mio1_27_w_512_1_1_512_d.dat
    $mycommand -np 512 ./main_1c_27_d.x  1 16 32 bench1w.inp > mio1_27_w_512_1_16_32_d.dat
    $mycommand -np 512 ./main_1c_27_d.x  4 8 16 bench1w.inp > mio1_27_w_512_4_8_16_d.dat
    $mycommand -np 512 ./main_1c_27_d.x  8 8 8 bench1w.inp > mio1_27_w_512_8_8_8_d.dat
    
    $mycommand -np 512 ./main_2c_27_d.x  1 1 512 benchw.inp > mio2_27_w_512_1_1_512_d.dat
    $mycommand -np 512 ./main_2c_27_d.x  1 16 32 benchw.inp > mio2_27_w_512_1_16_32_d.dat
    $mycommand -np 512 ./main_2c_27_d.x  4 8 16 benchw.inp > mio2_27_w_512_4_8_16_d.dat
    $mycommand -np 512 ./main_2c_27_d.x  8 8 8 benchw.inp > mio2_27_w_512_8_8_8_d.dat
    
    echo "Mode: STRONG - Running strong configuration"
    $mycommand -np 512 ./main_1c_27_d.x 8 8 8 bench1.inp > mio1_27_s_512_8_8_8_d.dat
    $mycommand -np 512 ./main_1c_27_d.x 4 8 16 bench1.inp > mio1_27_s_512_4_8_16_d.dat
    $mycommand -np 512 ./main_1c_27_d.x 1 16 32 bench1.inp > mio1_27_s_512_1_16_32_d.dat
    
    $mycommand -np 512 ./main_2c_27_d.x 8 8 8 bench.inp > mio2_27_s_512_8_8_8_d.dat
    $mycommand -np 512 ./main_2c_27_d.x 4 8 16 bench.inp > mio2_27_s_512_4_8_16_d.dat
    $mycommand -np 512 ./main_2c_27_d.x 1 16 32 bench.inp > mio2_27_s_512_1_16_32_d.dat

fi


if [[ "$mynodes" == "32" ]]; then
        echo "Mode: WEAK - Running weak configuration"

        $mycommand -np 128 ./main_1c_27high_d.x 1 1 128 bench1w.inp > mio1_27h_w_128_1_1_128_d.dat
        $mycommand -np 128 ./main_1c_27high_d.x 1 8 16 bench1w.inp > mio1_27h_w_128_1_8_16_d.dat
        $mycommand -np 128 ./main_1c_27high_d.x 4 4 8 bench1w.inp > mio1_27h_w_128_4_4_8_d.dat

	    $mycommand -np 128 ./main_2c_27high_d.x 1 1 128 benchw.inp > mio2_27h_w_128_1_1_128_d.dat
	    $mycommand -np 128 ./main_2c_27high_d.x 1 8 16 benchw.inp > mio2_27h_w_128_1_8_16_d.dat
        $mycommand -np 128 ./main_2c_27high_d.x 4 4 8 benchw.inp > mio2_27h_w_128_4_4_8_d.dat
	    
        echo "Mode: STRONG - Running strong configuration"
       $mycommand -np 128 ./main_1c_27high_d.x 4 4 8 bench1.inp > mio1_27h_s_128_4_4_8_d.dat
       $mycommand -np 128 ./main_1c_27high_d.x 1 8 16 bench1.inp > mio1_27h_s_128_1_8_16_d.dat

       $mycommand -np 128 ./main_2c_27high_d.x 4 4 8 bench.inp > mio2_27h_s_128_4_4_8_d.dat
       $mycommand -np 128 ./main_2c_27high_d.x 1 8 16 bench.inp > mio2_27h_s_128_1_8_16_d.dat

fi

if [[ "$mynodes" == "64" ]]; then

    echo "Mode: WEAK - Running weak configuration"
    
    $mycommand -np 256 ./main_1c_27high_d.x 1 1 256 bench1w.inp > mio1_27h_w_256_1_1_256_d.dat
    $mycommand -np 256 ./main_1c_27high_d.x 1 16 16 bench1w.inp > mio1_27h_w_256_1_16_16_d.dat
    $mycommand -np 256 ./main_1c_27high_d.x 4 8 8 bench1w.inp > mio1_27h_w_256_4_8_8_d.dat
    
    $mycommand -np 256 ./main_2c_27high_d.x 1 1 256 benchw.inp > mio2_27h_w_256_1_1_256_d.dat
    $mycommand -np 256 ./main_2c_27high_d.x 1 16 16 benchw.inp > mio2_27h_w_256_1_16_16_d.dat
    $mycommand -np 256 ./main_2c_27high_d.x 4 8 8 benchw.inp > mio2_27h_w_256_4_8_8_d.dat
    
    echo "Mode: STRONG - Running strong configuration"
    $mycommand -np 256 ./main_1c_27high_d.x 4 8 8 bench1.inp > mio1_27h_s_256_4_8_8_d.dat
    $mycommand -np 256 ./main_1c_27high_d.x 1 16 16 bench1.inp > mio1_27h_s_256_1_16_16_d.dat
    
    $mycommand -np 256 ./main_2c_27high_d.x 4 8 8 bench.inp > mio2_27h_s_256_4_8_8_d.dat
    $mycommand -np 256 ./main_2c_27high_d.x 1 16 16 bench.inp > mio2_27h_s_256_1_16_16_d.dat

fi

if [[ "$mynodes" == "128" ]]; then
    echo "Mode: WEAK - Running weak configuration"

    $mycommand -np 512 ./main_1c_27high_d.x  1 1 512 bench1w.inp > mio1_27h_w_512_1_1_512_d.dat
    $mycommand -np 512 ./main_1c_27high_d.x  1 16 32 bench1w.inp > mio1_27h_w_512_1_16_32_d.dat
    $mycommand -np 512 ./main_1c_27high_d.x  4 8 16 bench1w.inp > mio1_27h_w_512_4_8_16_d.dat
    $mycommand -np 512 ./main_1c_27high_d.x  8 8 8 bench1w.inp > mio1_27h_w_512_8_8_8_d.dat
    
    $mycommand -np 512 ./main_2c_27high_d.x  1 1 512 benchw.inp > mio2_27h_w_512_1_1_512_d.dat
    $mycommand -np 512 ./main_2c_27high_d.x  1 16 32 benchw.inp > mio2_27h_w_512_1_16_32_d.dat
    $mycommand -np 512 ./main_2c_27high_d.x  4 8 16 benchw.inp > mio2_27h_w_512_4_8_16_d.dat
    $mycommand -np 512 ./main_2c_27high_d.x  8 8 8 benchw.inp > mio2_27h_w_512_8_8_8_d.dat
    
    echo "Mode: STRONG - Running strong configuration"
    $mycommand -np 512 ./main_1c_27high_d.x 8 8 8 bench1.inp > mio1_27h_s_512_8_8_8_d.dat
    $mycommand -np 512 ./main_1c_27high_d.x 4 8 16 bench1.inp > mio1_27h_s_512_4_8_16_d.dat
    $mycommand -np 512 ./main_1c_27high_d.x 1 16 32 bench1.inp > mio1_27h_s_512_1_16_32_d.dat
    
    $mycommand -np 512 ./main_2c_27high_d.x 8 8 8 bench.inp > mio2_27h_s_512_8_8_8_d.dat
    $mycommand -np 512 ./main_2c_27high_d.x 4 8 16 bench.inp > mio2_27h_s_512_4_8_16_d.dat
    $mycommand -np 512 ./main_2c_27high_d.x 1 16 32 bench.inp > mio2_27h_s_512_1_16_32_d.dat

fi



##########################################################################################################################



if [[ "$mynodes" == "32" ]]; then
        echo "Mode: WEAK - Running weak configuration"

        $mycommand -np 128 ./main_1c_15_sd.x 1 1 128 bench1w.inp > mio1_15_w_128_1_1_128_sd.dat
        $mycommand -np 128 ./main_1c_15_sd.x 1 8 16 bench1w.inp > mio1_15_w_128_1_8_16_sd.dat
        $mycommand -np 128 ./main_1c_15_sd.x 4 4 8 bench1w.inp > mio1_15_w_128_4_4_8_sd.dat

	    $mycommand -np 128 ./main_2c_15_sd.x 1 1 128 benchw.inp > mio2_15_w_128_1_1_128_sd.dat
	    $mycommand -np 128 ./main_2c_15_sd.x 1 8 16 benchw.inp > mio2_15_w_128_1_8_16_sd.dat
        $mycommand -np 128 ./main_2c_15_sd.x 4 4 8 benchw.inp > mio2_15_w_128_4_4_8_sd.dat
	    
        echo "Mode: STRONG - Running strong configuration"
       $mycommand -np 128 ./main_1c_15_sd.x 4 4 8 bench1.inp > mio1_15_s_128_4_4_8_sd.dat
       $mycommand -np 128 ./main_1c_15_sd.x 1 8 16 bench1.inp > mio1_15_s_128_1_8_16_sd.dat

       $mycommand -np 128 ./main_2c_15_sd.x 4 4 8 bench.inp > mio2_15_s_128_4_4_8_sd.dat
       $mycommand -np 128 ./main_2c_15_sd.x 1 8 16 bench.inp > mio2_15_s_128_1_8_16_sd.dat

fi

if [[ "$mynodes" == "64" ]]; then

    echo "Mode: WEAK - Running weak configuration"
    
    $mycommand -np 256 ./main_1c_15_sd.x 1 1 256 bench1w.inp > mio1_15_w_256_1_1_256_sd.dat
    $mycommand -np 256 ./main_1c_15_sd.x 1 16 16 bench1w.inp > mio1_15_w_256_1_16_16_sd.dat
    $mycommand -np 256 ./main_1c_15_sd.x 4 8 8 bench1w.inp > mio1_15_w_256_4_8_8_sd.dat
    
    $mycommand -np 256 ./main_2c_15_sd.x 1 1 256 benchw.inp > mio2_15_w_256_1_1_256_sd.dat
    $mycommand -np 256 ./main_2c_15_sd.x 1 16 16 benchw.inp > mio2_15_w_256_1_16_16_sd.dat
    $mycommand -np 256 ./main_2c_15_sd.x 4 8 8 benchw.inp > mio2_15_w_256_4_8_8_sd.dat
    
    echo "Mode: STRONG - Running strong configuration"
    $mycommand -np 256 ./main_1c_15_sd.x 4 8 8 bench1.inp > mio1_15_s_256_4_8_8_sd.dat
    $mycommand -np 256 ./main_1c_15_sd.x 1 16 16 bench1.inp > mio1_15_s_256_1_16_16_sd.dat
    
    $mycommand -np 256 ./main_2c_15_sd.x 4 8 8 bench.inp > mio2_15_s_256_4_8_8_sd.dat
    $mycommand -np 256 ./main_2c_15_sd.x 1 16 16 bench.inp > mio2_15_s_256_1_16_16_sd.dat

fi

if [[ "$mynodes" == "128" ]]; then
    echo "Mode: WEAK - Running weak configuration"
    $mycommand -np 512 ./main_1c_15_sd.x  1 1 512 bench1w.inp > mio1_15_w_512_1_1_512_sd.dat
    $mycommand -np 512 ./main_1c_15_sd.x  1 16 32 bench1w.inp > mio1_15_w_512_1_16_32_sd.dat
    $mycommand -np 512 ./main_1c_15_sd.x  4 8 16 bench1w.inp > mio1_15_w_512_4_8_16_sd.dat
    $mycommand -np 512 ./main_1c_15_sd.x  8 8 8 bench1w.inp > mio1_15_w_512_8_8_8_sd.dat
    
    $mycommand -np 512 ./main_2c_15_sd.x  1 1 512 benchw.inp > mio2_15_w_512_1_1_512_sd.dat
    $mycommand -np 512 ./main_2c_15_sd.x  1 16 32 benchw.inp > mio2_15_w_512_1_16_32_sd.dat
    $mycommand -np 512 ./main_2c_15_sd.x  4 8 16 benchw.inp > mio2_15_w_512_4_8_16_sd.dat
    $mycommand -np 512 ./main_2c_15_sd.x  8 8 8 benchw.inp > mio2_15_w_512_8_8_8_sd.dat
    
    echo "Mode: STRONG - Running strong configuration"
    $mycommand -np 512 ./main_1c_15_sd.x 8 8 8 bench1.inp > mio1_15_s_512_8_8_8_sd.dat
    $mycommand -np 512 ./main_1c_15_sd.x 4 8 16 bench1.inp > mio1_15_s_512_4_8_16_sd.dat
    $mycommand -np 512 ./main_1c_15_sd.x 1 16 32 bench1.inp > mio1_15_s_512_1_16_32_sd.dat
    
    $mycommand -np 512 ./main_2c_15_sd.x 8 8 8 bench.inp > mio2_15_s_512_8_8_8_sd.dat
    $mycommand -np 512 ./main_2c_15_sd.x 4 8 16 bench.inp > mio2_15_s_512_4_8_16_sd.dat
    $mycommand -np 512 ./main_2c_15_sd.x 1 16 32 bench.inp > mio2_15_s_512_1_16_32_sd.dat

fi

if [[ "$mynodes" == "32" ]]; then
        echo "Mode: WEAK - Running weak configuration"
       
        $mycommand -np 128 ./main_1c_19_sd.x 1 1 128 bench1w.inp > mio1_19_w_128_1_1_128_sd.dat
        $mycommand -np 128 ./main_1c_19_sd.x 1 8 16 bench1w.inp > mio1_19_w_128_1_8_16_sd.dat
        $mycommand -np 128 ./main_1c_19_sd.x 4 4 8 bench1w.inp > mio1_19_w_128_4_4_8_sd.dat

	    $mycommand -np 128 ./main_2c_19_sd.x 1 1 128 benchw.inp > mio2_19_w_128_1_1_128_sd.dat
	    $mycommand -np 128 ./main_2c_19_sd.x 1 8 16 benchw.inp > mio2_19_w_128_1_8_16_sd.dat
        $mycommand -np 128 ./main_2c_19_sd.x 4 4 8 benchw.inp > mio2_19_w_128_4_4_8_sd.dat
	    
        echo "Mode: STRONG - Running strong configuration"
       $mycommand -np 128 ./main_1c_19_sd.x 4 4 8 bench1.inp > mio1_19_s_128_4_4_8_sd.dat
       $mycommand -np 128 ./main_1c_19_sd.x 1 8 16 bench1.inp > mio1_19_s_128_1_8_16_sd.dat

       $mycommand -np 128 ./main_2c_19_sd.x 4 4 8 bench.inp > mio2_19_s_128_4_4_8_sd.dat
       $mycommand -np 128 ./main_2c_19_sd.x 1 8 16 bench.inp > mio2_19_s_128_1_8_16_sd.dat

fi

if [[ "$mynodes" == "64" ]]; then

    echo "Mode: WEAK - Running weak configuration"
    
    $mycommand -np 256 ./main_1c_19_sd.x 1 1 256 bench1w.inp > mio1_19_w_256_1_1_256_sd.dat
    $mycommand -np 256 ./main_1c_19_sd.x 1 16 16 bench1w.inp > mio1_19_w_256_1_16_16_sd.dat
    $mycommand -np 256 ./main_1c_19_sd.x 4 8 8 bench1w.inp > mio1_19_w_256_4_8_8_sd.dat
    
    $mycommand -np 256 ./main_2c_19_sd.x 1 1 256 benchw.inp > mio2_19_w_256_1_1_256_sd.dat
    $mycommand -np 256 ./main_2c_19_sd.x 1 16 16 benchw.inp > mio2_19_w_256_1_16_16_sd.dat
    $mycommand -np 256 ./main_2c_19_sd.x 4 8 8 benchw.inp > mio2_19_w_256_4_8_8_sd.dat
    
    echo "Mode: STRONG - Running strong configuration"
    $mycommand -np 256 ./main_1c_19_sd.x 4 8 8 bench1.inp > mio1_19_s_256_4_8_8_sd.dat
    $mycommand -np 256 ./main_1c_19_sd.x 1 16 16 bench1.inp > mio1_19_s_256_1_16_16_sd.dat
    
    $mycommand -np 256 ./main_2c_19_sd.x 4 8 8 bench.inp > mio2_19_s_256_4_8_8_sd.dat
    $mycommand -np 256 ./main_2c_19_sd.x 1 16 16 bench.inp > mio2_19_s_256_1_16_16_sd.dat

fi

if [[ "$mynodes" == "128" ]]; then
    echo "Mode: WEAK - Running weak configuration"
    
    $mycommand -np 512 ./main_1c_19_sd.x  1 1 512 bench1w.inp > mio1_19_w_512_1_1_512_sd.dat
    $mycommand -np 512 ./main_1c_19_sd.x  1 16 32 bench1w.inp > mio1_19_w_512_1_16_32_sd.dat
    $mycommand -np 512 ./main_1c_19_sd.x  4 8 16 bench1w.inp > mio1_19_w_512_4_8_16_sd.dat
    $mycommand -np 512 ./main_1c_19_sd.x  8 8 8 bench1w.inp > mio1_19_w_512_8_8_8_sd.dat
    
    $mycommand -np 512 ./main_2c_19_sd.x  1 1 512 benchw.inp > mio2_19_w_512_1_1_512_sd.dat
    $mycommand -np 512 ./main_2c_19_sd.x  1 16 32 benchw.inp > mio2_19_w_512_1_16_32_sd.dat
    $mycommand -np 512 ./main_2c_19_sd.x  4 8 16 benchw.inp > mio2_19_w_512_4_8_16_sd.dat
    $mycommand -np 512 ./main_2c_19_sd.x  8 8 8 benchw.inp > mio2_19_w_512_8_8_8_sd.dat
    
    echo "Mode: STRONG - Running strong configuration"
    $mycommand -np 512 ./main_1c_19_sd.x 8 8 8 bench1.inp > mio1_19_s_512_8_8_8_sd.dat
    $mycommand -np 512 ./main_1c_19_sd.x 4 8 16 bench1.inp > mio1_19_s_512_4_8_16_sd.dat
    $mycommand -np 512 ./main_1c_19_sd.x 1 16 32 bench1.inp > mio1_19_s_512_1_16_32_sd.dat
    
    $mycommand -np 512 ./main_2c_19_sd.x 8 8 8 bench.inp > mio2_19_s_512_8_8_8_sd.dat
    $mycommand -np 512 ./main_2c_19_sd.x 4 8 16 bench.inp > mio2_19_s_512_4_8_16_sd.dat
    $mycommand -np 512 ./main_2c_19_sd.x 1 16 32 bench.inp > mio2_19_s_512_1_16_32_sd.dat


fi


if [[ "$mynodes" == "32" ]]; then
        echo "Mode: WEAK - Running weak configuration"

        $mycommand -np 128 ./main_1c_27_sd.x 1 1 128 bench1w.inp > mio1_27_w_128_1_1_128_sd.dat
        $mycommand -np 128 ./main_1c_27_sd.x 1 8 16 bench1w.inp > mio1_27_w_128_1_8_16_sd.dat
        $mycommand -np 128 ./main_1c_27_sd.x 4 4 8 bench1w.inp > mio1_27_w_128_4_4_8_sd.dat

	    $mycommand -np 128 ./main_2c_27_sd.x 1 1 128 benchw.inp > mio2_27_w_128_1_1_128_sd.dat
	    $mycommand -np 128 ./main_2c_27_sd.x 1 8 16 benchw.inp > mio2_27_w_128_1_8_16_sd.dat
        $mycommand -np 128 ./main_2c_27_sd.x 4 4 8 benchw.inp > mio2_27_w_128_4_4_8_sd.dat
	    
        echo "Mode: STRONG - Running strong configuration"
       $mycommand -np 128 ./main_1c_27_sd.x 4 4 8 bench1.inp > mio1_27_s_128_4_4_8_sd.dat
       $mycommand -np 128 ./main_1c_27_sd.x 1 8 16 bench1.inp > mio1_27_s_128_1_8_16_sd.dat

       $mycommand -np 128 ./main_2c_27_sd.x 4 4 8 bench.inp > mio2_27_s_128_4_4_8_sd.dat
       $mycommand -np 128 ./main_2c_27_sd.x 1 8 16 bench.inp > mio2_27_s_128_1_8_16_sd.dat

fi

if [[ "$mynodes" == "64" ]]; then

    echo "Mode: WEAK - Running weak configuration"
    
    $mycommand -np 256 ./main_1c_27_sd.x 1 1 256 bench1w.inp > mio1_27_w_256_1_1_256_sd.dat
    $mycommand -np 256 ./main_1c_27_sd.x 1 16 16 bench1w.inp > mio1_27_w_256_1_16_16_sd.dat
    $mycommand -np 256 ./main_1c_27_sd.x 4 8 8 bench1w.inp > mio1_27_w_256_4_8_8_sd.dat
    
    $mycommand -np 256 ./main_2c_27_sd.x 1 1 256 benchw.inp > mio2_27_w_256_1_1_256_sd.dat
    $mycommand -np 256 ./main_2c_27_sd.x 1 16 16 benchw.inp > mio2_27_w_256_1_16_16_sd.dat
    $mycommand -np 256 ./main_2c_27_sd.x 4 8 8 benchw.inp > mio2_27_w_256_4_8_8_sd.dat
    
    echo "Mode: STRONG - Running strong configuration"
    $mycommand -np 256 ./main_1c_27_sd.x 4 8 8 bench1.inp > mio1_27_s_256_4_8_8_sd.dat
    $mycommand -np 256 ./main_1c_27_sd.x 1 16 16 bench1.inp > mio1_27_s_256_1_16_16_sd.dat
    
    $mycommand -np 256 ./main_2c_27_sd.x 4 8 8 bench.inp > mio2_27_s_256_4_8_8_sd.dat
    $mycommand -np 256 ./main_2c_27_sd.x 1 16 16 bench.inp > mio2_27_s_256_1_16_16_sd.dat

fi

if [[ "$mynodes" == "128" ]]; then
    echo "Mode: WEAK - Running weak configuration"

    $mycommand -np 512 ./main_1c_27_sd.x  1 1 512 bench1w.inp > mio1_27_w_512_1_1_512_sd.dat
    $mycommand -np 512 ./main_1c_27_sd.x  1 16 32 bench1w.inp > mio1_27_w_512_1_16_32_sd.dat
    $mycommand -np 512 ./main_1c_27_sd.x  4 8 16 bench1w.inp > mio1_27_w_512_4_8_16_sd.dat
    $mycommand -np 512 ./main_1c_27_sd.x  8 8 8 bench1w.inp > mio1_27_w_512_8_8_8_sd.dat
    
    $mycommand -np 512 ./main_2c_27_sd.x  1 1 512 benchw.inp > mio2_27_w_512_1_1_512_sd.dat
    $mycommand -np 512 ./main_2c_27_sd.x  1 16 32 benchw.inp > mio2_27_w_512_1_16_32_sd.dat
    $mycommand -np 512 ./main_2c_27_sd.x  4 8 16 benchw.inp > mio2_27_w_512_4_8_16_sd.dat
    $mycommand -np 512 ./main_2c_27_sd.x  8 8 8 benchw.inp > mio2_27_w_512_8_8_8_sd.dat
    
    echo "Mode: STRONG - Running strong configuration"
    $mycommand -np 512 ./main_1c_27_sd.x 8 8 8 bench1.inp > mio1_27_s_512_8_8_8_sd.dat
    $mycommand -np 512 ./main_1c_27_sd.x 4 8 16 bench1.inp > mio1_27_s_512_4_8_16_sd.dat
    $mycommand -np 512 ./main_1c_27_sd.x 1 16 32 bench1.inp > mio1_27_s_512_1_16_32_sd.dat
    
    $mycommand -np 512 ./main_2c_27_sd.x 8 8 8 bench.inp > mio2_27_s_512_8_8_8_sd.dat
    $mycommand -np 512 ./main_2c_27_sd.x 4 8 16 bench.inp > mio2_27_s_512_4_8_16_sd.dat
    $mycommand -np 512 ./main_2c_27_sd.x 1 16 32 bench.inp > mio2_27_s_512_1_16_32_sd.dat

fi


if [[ "$mynodes" == "32" ]]; then
        echo "Mode: WEAK - Running weak configuration"

        $mycommand -np 128 ./main_1c_27high_sd.x 1 1 128 bench1w.inp > mio1_27h_w_128_1_1_128_sd.dat
        $mycommand -np 128 ./main_1c_27high_sd.x 1 8 16 bench1w.inp > mio1_27h_w_128_1_8_16_sd.dat
        $mycommand -np 128 ./main_1c_27high_sd.x 4 4 8 bench1w.inp > mio1_27h_w_128_4_4_8_sd.dat

	    $mycommand -np 128 ./main_2c_27high_sd.x 1 1 128 benchw.inp > mio2_27h_w_128_1_1_128_sd.dat
	    $mycommand -np 128 ./main_2c_27high_sd.x 1 8 16 benchw.inp > mio2_27h_w_128_1_8_16_sd.dat
        $mycommand -np 128 ./main_2c_27high_sd.x 4 4 8 benchw.inp > mio2_27h_w_128_4_4_8_sd.dat
	    
        echo "Mode: STRONG - Running strong configuration"
       $mycommand -np 128 ./main_1c_27high_sd.x 4 4 8 bench1.inp > mio1_27h_s_128_4_4_8_sd.dat
       $mycommand -np 128 ./main_1c_27high_sd.x 1 8 16 bench1.inp > mio1_27h_s_128_1_8_16_sd.dat

       $mycommand -np 128 ./main_2c_27high_sd.x 4 4 8 bench.inp > mio2_27h_s_128_4_4_8_sd.dat
       $mycommand -np 128 ./main_2c_27high_sd.x 1 8 16 bench.inp > mio2_27h_s_128_1_8_16_sd.dat

fi

if [[ "$mynodes" == "64" ]]; then

    echo "Mode: WEAK - Running weak configuration"
    
    $mycommand -np 256 ./main_1c_27high_sd.x 1 1 256 bench1w.inp > mio1_27h_w_256_1_1_256_sd.dat
    $mycommand -np 256 ./main_1c_27high_sd.x 1 16 16 bench1w.inp > mio1_27h_w_256_1_16_16_sd.dat
    $mycommand -np 256 ./main_1c_27high_sd.x 4 8 8 bench1w.inp > mio1_27h_w_256_4_8_8_sd.dat
    
    $mycommand -np 256 ./main_2c_27high_sd.x 1 1 256 benchw.inp > mio2_27h_w_256_1_1_256_sd.dat
    $mycommand -np 256 ./main_2c_27high_sd.x 1 16 16 benchw.inp > mio2_27h_w_256_1_16_16_sd.dat
    $mycommand -np 256 ./main_2c_27high_sd.x 4 8 8 benchw.inp > mio2_27h_w_256_4_8_8_sd.dat
    
    echo "Mode: STRONG - Running strong configuration"
    $mycommand -np 256 ./main_1c_27high_sd.x 4 8 8 bench1.inp > mio1_27h_s_256_4_8_8_sd.dat
    $mycommand -np 256 ./main_1c_27high_sd.x 1 16 16 bench1.inp > mio1_27h_s_256_1_16_16_sd.dat
    
    $mycommand -np 256 ./main_2c_27high_sd.x 4 8 8 bench.inp > mio2_27h_s_256_4_8_8_sd.dat
    $mycommand -np 256 ./main_2c_27high_sd.x 1 16 16 bench.inp > mio2_27h_s_256_1_16_16_sd.dat

fi

if [[ "$mynodes" == "128" ]]; then
    echo "Mode: WEAK - Running weak configuration"

    $mycommand -np 512 ./main_1c_27high_sd.x  1 1 512 bench1w.inp > mio1_27h_w_512_1_1_512_sd.dat
    $mycommand -np 512 ./main_1c_27high_sd.x  1 16 32 bench1w.inp > mio1_27h_w_512_1_16_32_sd.dat
    $mycommand -np 512 ./main_1c_27high_sd.x  4 8 16 bench1w.inp > mio1_27h_w_512_4_8_16_sd.dat
    $mycommand -np 512 ./main_1c_27high_sd.x  8 8 8 bench1w.inp > mio1_27h_w_512_8_8_8_sd.dat
    
    $mycommand -np 512 ./main_2c_27high_sd.x  1 1 512 benchw.inp > mio2_27h_w_512_1_1_512_sd.dat
    $mycommand -np 512 ./main_2c_27high_sd.x  1 16 32 benchw.inp > mio2_27h_w_512_1_16_32_sd.dat
    $mycommand -np 512 ./main_2c_27high_sd.x  4 8 16 benchw.inp > mio2_27h_w_512_4_8_16_sd.dat
    $mycommand -np 512 ./main_2c_27high_sd.x  8 8 8 benchw.inp > mio2_27h_w_512_8_8_8_sd.dat
    
    echo "Mode: STRONG - Running strong configuration"
    $mycommand -np 512 ./main_1c_27high_sd.x 8 8 8 bench1.inp > mio1_27h_s_512_8_8_8_sd.dat
    $mycommand -np 512 ./main_1c_27high_sd.x 4 8 16 bench1.inp > mio1_27h_s_512_4_8_16_sd.dat
    $mycommand -np 512 ./main_1c_27high_sd.x 1 16 32 bench1.inp > mio1_27h_s_512_1_16_32_sd.dat
    
    $mycommand -np 512 ./main_2c_27high_sd.x 8 8 8 bench.inp > mio2_27h_s_512_8_8_8_sd.dat
    $mycommand -np 512 ./main_2c_27high_sd.x 4 8 16 bench.inp > mio2_27h_s_512_4_8_16_sd.dat
    $mycommand -np 512 ./main_2c_27high_sd.x 1 16 32 bench.inp > mio2_27h_s_512_1_16_32_sd.dat

fi

