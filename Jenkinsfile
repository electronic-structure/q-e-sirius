#!groovy

pipeline {
    agent {
        node {
            label 'ssl_daintvm1'
        }
    }
    //environment {
    //    EB_CUSTOM_REPOSITORY = '/users/simonpi/jenkins/production/easybuild'
    //}
    stages {
        stage('Checkout') {
            steps {
                dir('q-e-sirius') {
                    checkout scm
                    echo "Running ${env.BUILD_ID} on ${env.JENKINS_URL}"
                }
            }
        }
        //stage('Compile') {
        //    steps {
        //        dir('SIRIUS') {
        //            sh '''
        //                   #!/bin/bash -l
        //                   export ENVFILE=$(realpath ci/env-gnu-gpu)
        //                   rm -f build-daint-gpu.out
        //                   rm -f build-daint-gpu.err
        //                   sbatch --wait ci/build-daint-gpu.sh
        //                   echo "---------- build-daint-gpu.out ----------"
        //                   cat build-daint-gpu.out
        //                   echo "---------- build-daint-gpu.err ----------"
        //                   cat build-daint-gpu.err
        //                   # check that sirius.scf has been built
        //                   type -f build/apps/dft_loop/sirius.scf
        //                   '''
        //        }
        //    }
        //}
    }
    post {
        always {
            archiveArtifacts artifacts: '**/*.out', fingerprint: true
            archiveArtifacts artifacts: '**/*.err', fingerprint: true
        }
    }
}

