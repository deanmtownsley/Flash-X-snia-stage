node{
    
    stage('Organize'){
        sh "rm -rf ./*"
	
        
        dir('FLASH-X'){
            checkout scm
            sh 'cp /home/jenkins/Makefile.h_flash5 sites/Prototypes/Linux/Makefile.h'
            sh 'perl -i -p -e \'s|git@(.*?):(.*?)|https://\\1/\\2|g\' .gitmodules'
            sh 'git submodule update --init'
            sh 'ls'
	    dir('tools/sfocu'){
	        sh 'cp Makefile.ganon Makefile'
		sh 'SITE=Prototypes/Linux make'
	    }
        }
    }

    stage('DeleptonizationWave'){ catchError(buildResult: 'UNSTABLE', stageResult: 'FAILURE'){
	dir('FLASH-X'){
	    sh './setup DeleptonizationWave -auto -1d +spherical +mode1 -nxb=16 -maxblocks=100 nE=5 nSpecies=2 nNodes=2 nMoments=4 momentClosure=MINERBO thornadoSolver=EMAB thornadoOrder=ORDER_1 -with-unit=source/physics/Eos/EosMain/WeakLib +uhd +pm4dev -parfile=test_sph_1d.par -objdir=dele1d_MI_5E'
	    dir('dele1d_MI_5E'){
	        sh 'make -j'
	    	sh 'ln -s /home/jenkins/wl-Op-SFHo-15-25-50-E40-B85-AbEm.h5'
		sh 'ln -s /home/jenkins/wl-Op-SFHo-15-25-50-E40-B85-Iso.h5'
		sh 'ln -s /home/jenkins/wl-Op-SFHo-15-25-50-E40-B85-NES.h5'
		sh 'ln -s /home/jenkins/wl-Op-SFHo-15-25-50-E40-B85-Pair.h5'
		sh './flash5'
		sh '../tools/sfocu/sfocu /home/jenkins/Sep06_NewDelep_baseline/deleptonizationwave_hdf5_chk_0001 deleptonizationwave_hdf5_chk_0001 > DeleptonizationWave.txt || echo "Return Code: $?" '
		archiveArtifacts: 'DeleptonizationWave.txt'
	    }
	}
    }}

    stage('CartRadEqm-photon'){catchError(buildResult: 'UNSTABLE', stageResult: 'FAILURE'){
        dir('FLASH-X'){
            sh 'python3 bin/setup.py CartRadEqm -auto -3d -maxblocks=200 -debug -objdir=cartradeqm-photon Bittree=True FlashAvoidOrrery=True -with-unit=source/Particles/ParticlesMain/active/ParticlesOwned/MonteCarlo/PhotonRad'
            dir('cartradeqm-photon'){
                sh 'make -j'
                sh 'mpirun -np 8 ./flash5'
                sh 'python3 plot_temperature.py'
                archiveArtifacts artifacts: '*.pdf'
            }
        }
    }}

    
    stage('CartRadDiff-photon'){catchError(buildResult: 'UNSTABLE', stageResult: 'FAILURE'){
        dir('FLASH-X'){
            sh 'python3 bin/setup.py CartRadDiff -auto -3d -maxblocks=200 -debug -objdir=cartraddiff-photon Bittree=True FlashAvoidOrrery=True -with-unit=source/Particles/ParticlesMain/active/ParticlesOwned/MonteCarlo/PhotonRad'
            dir('cartraddiff-photon'){
                sh 'make -j'
                sh 'mpirun -np 8 ./flash5'
                sh 'python3 urad_from_mcps.py cartimc_1M_hdf5_chk_ 1 3'
                archiveArtifacts artifacts: '*.pdf'
            }
        }
    }}

    stage('CartRadEqm-neutrino'){catchError(buildResult: 'UNSTABLE', stageResult: 'FAILURE'){
        dir('FLASH-X'){
            sh 'python3 bin/setup.py CartRadEqm -auto -3d -nxb=2 -nyb=2 -nzb=2 -maxblocks=200 nE=20 nSpecies=2 -debug -objdir=cartradeqm-neutrino Bittree=True FlashAvoidOrrery=True -with-unit=source/physics/Eos/EosMain/WeakLib -with-unit=source/Particles/ParticlesMain/active/ParticlesOwned/MonteCarlo/NeutrinoRad'
            dir('cartradeqm-neutrino'){
                sh 'make -j'
		sh 'ln -s /home/jenkins/wl-Op-SFHo-15-25-50-E40-B85-AbEm.h5'
		sh 'ln -s /home/jenkins/wl-Op-SFHo-15-25-50-E40-B85-Iso.h5'
		sh 'ln -s /home/jenkins/wl-Op-SFHo-15-25-50-E40-B85-NES.h5'
		sh 'ln -s /home/jenkins/wl-Op-SFHo-15-25-50-E40-B85-Pair.h5'
		sh 'cp ../source/Simulation/SimulationMain/CartRadEqm/flash.par_neutrino flash.par'
                sh 'mpirun -np 8 ./flash5'
                sh 'python3 plot_temperature.py'
                archiveArtifacts artifacts: '*.pdf'
            }
        }
    }}

    stage('CartRadDiff-neutrino'){catchError(buildResult: 'UNSTABLE', stageResult: 'FAILURE'){
        dir('FLASH-X'){
            sh 'python3 bin/setup.py CartRadDiff -auto -3d nE=20 nSpecies=2 -maxblocks=200 -debug -objdir=cartraddiff-neutrino Bittree=True FlashAvoidOrrery=True -with-unit=source/physics/Eos/EosMain/WeakLib -with-unit=source/Particles/ParticlesMain/active/ParticlesOwned/MonteCarlo/NeutrinoRad'
            dir('cartraddiff-neutrino'){
                sh 'make -j'
		sh 'ln -s /home/jenkins/wl-Op-SFHo-15-25-50-E40-B85-AbEm.h5'
		sh 'ln -s /home/jenkins/wl-Op-SFHo-15-25-50-E40-B85-Iso.h5'
		sh 'ln -s /home/jenkins/wl-Op-SFHo-15-25-50-E40-B85-NES.h5'
		sh 'ln -s /home/jenkins/wl-Op-SFHo-15-25-50-E40-B85-Pair.h5'
		sh 'cp ../source/Simulation/SimulationMain/CartRadDiff/flash.par_neutrinos flash.par'
                sh 'mpirun -np 8 ./flash5'
                sh 'python3 urad_from_mcps.py cartimc_1M_hdf5_chk_ 1 3'
                archiveArtifacts artifacts: '*.pdf'
            }
        }
    }}
    
    stage('RadiativeShock-photon'){catchError(buildResult: 'UNSTABLE', stageResult: 'FAILURE'){
        dir('FLASH-X'){
            sh 'python3 bin/setup.py RadiativeShock -auto -1d -maxblocks=100 -objdir=radiativeshock1d-photon -debug Bittree=True FlashAvoidOrrery=True -with-unit=source/Particles/ParticlesMain/active/ParticlesOwned/MonteCarlo/PhotonRad'
            dir('radiativeshock1d-photon'){
                sh 'make -j'
                sh 'mpirun -np 8 ./flash5'
		sh 'python3 radshock_super_1d.py'
		archiveArtifacts artifacts: '*.pdf'
            }
        }
    }}
}
