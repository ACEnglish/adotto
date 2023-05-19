sbatch -J refine_hipstr -o logs/refine_hipstr.log -e logs/refine_hipstr.log --time 500 --partition medium --nodes=1 --ntasks-per-node=16 --mem=24000 jobs/refine_hipstr.sh
sbatch -J refine_gangstr -o logs/refine_gangstr.log -e logs/refine_gangstr.log --time 500 --partition medium --nodes=1 --ntasks-per-node=16 --mem=24000 jobs/refine_gangstr.sh
sbatch -J refine_trgt -o logs/refine_trgt.log -e logs/refine_trgt.log --time 500 --partition medium --nodes=1 --ntasks-per-node=16 --mem=24000 jobs/refine_trgt.sh
sbatch -J refine_biograph -o logs/refine_biograph.log -e logs/refine_biograph.log --time 500 --partition medium --nodes=1 --ntasks-per-node=16 --mem=24000 jobs/refine_biograph.sh
