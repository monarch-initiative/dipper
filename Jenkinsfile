/**
    Consider linting before running
    see:  https://jenkins.io/doc/book/pipeline/development/
    curl -X POST -H $(curl "127.0.0.1/crumbIssuer/api/xml?xpath=concat(//crumbRequestField,\":\",//crumb)") -F "jenkinsfile=<Jenkinsfile" 127.0.0.1/pipeline-model-converter/validate

**/

pipeline {

    agent any

    /*triggers {
         Run every Monday at 5pm
         cron('H 17 * * 1')
    }*/

    environment {
        // Pin dates and day to beginning of run.
        START_DATE = sh(
            script: 'date +%Y-%m-%d',
            returnStdout: true
        ).trim()

        START_DAY = sh(
            script: 'date +%A',
            returnStdout: true
        ).trim()

        YYYYMM = sh(
            script: 'date +%Y%m',
            returnStdout: true
        ).trim()

        MONARCHIVE = 'monarch@monarch-archive:/var/www/data/$YYYYMM/'
        DIPPERCACHE = 'https://archive.monarchinitiative.org/DipperCache'
        MONARCH_DATA_FS = 'monarch@monarch-ttl-prod'
        DIPPER = "venv/bin/python dipper-etl.py --skip_tests --data_release_version $YYYYMM"

        // https://issues.jenkins-ci.org/browse/JENKINS-47881
        DATA_DEST = "${env.RELEASE ? '/var/www/data/dev/' : '/var/www/data/experimental/'}"
        MONARCH_DATA_DEST = "$MONARCH_DATA_FS:$DATA_DEST"

        /* human, mouse, zebrafish, fly, worm, frog */
        COMMON_TAXON = "9606,10090,7955,7227,6239,8364"
        /* 10116 is rat and might be included if found relevant where it is now missing */

    }

    options {
        buildDiscarder(logRotator(numToKeepStr: '14'))
    }

    stages {
        stage('Build dipper package') {
            steps {
                dir('./config') {
                    git(
                        url: 'https://github.com/monarch-initiative/configs.git',
                        credentialsId: 'c36e48af-8e33-4977-be25-0555ebb4275a',
                        branch: 'master'
                    )
                    sh '''
                        cd .. && cp config/Dipper/conf.yaml ./dipper/conf.yaml
                        virtualenv -p /usr/bin/python3.6 venv
                        venv/bin/pip install -r requirements.txt
                        venv/bin/pip install -r requirements/all-sources.txt

                        echo "Clean up previous runs"
                        rm -f ./out/*.ttl ./out/*.nt
                        echo "Anything remaining should not still be in './out'"
                        rm -fr ./out
                        rm -fr ./raw
                    '''
                }
            }
        }
        stage("Validate Jenkinsfile"){
            steps{
                sh '''
                    curl -s -X POST -H $(curl "127.0.0.1/crumbIssuer/api/xml?xpath=concat(//crumbRequestField,\":\",//crumb)") -F "jenkinsfile=<Jenkinsfile" 127.0.0.1/pipeline-model-converter/validate
                '''
            }
        }
        stage('Generate monarch rdf') {
            parallel {
                stage("Process sources that call OMIM") {
                    stages {
                        stage("OMIM") {
                            when {
                                anyOf {
                                    expression { env.RUN_ALL != null }
                                    expression { env.OMIM != null }
                                }
                            }
                            steps {
                                sh '''
                                    SOURCE=omim
                                    $DIPPER --sources $SOURCE --quiet
                                    scp ./out/${SOURCE}.ttl ./out/${SOURCE}_dataset.ttl $MONARCH_DATA_DEST
                                '''
                            }
                        }
                        stage("NCBI Gene") {
                            when {
                                anyOf {
                                    expression { env.RUN_ALL != null }
                                    expression { env.NCBIGENE != null }
                                }
                            }
                            steps {
                                sh '''
                                    SOURCE=ncbigene
                                    $DIPPER --sources $SOURCE \
                                        --taxon $COMMON_TAXON,10116,28377,3702,9913,9615,9031,44689,9796,9544,13616,9258,9598,9823,4896,31033,8364,9685,559292,8355
                                    scp ./out/${SOURCE}.ttl ./out/${SOURCE}_dataset.ttl $MONARCH_DATA_DEST
                                '''
                            }
                        }
                        stage("HGNC") {
                            when {
                                anyOf {
                                    expression { env.RUN_ALL != null }
                                    expression { env.HGNC != null }
                                }
                            }
                            steps {
                                sh '''
                                    SOURCE=hgnc
                                    $DIPPER --sources $SOURCE
                                    scp ./out/${SOURCE}.ttl ./out/${SOURCE}_dataset.ttl $MONARCH_DATA_DEST
                                '''
                            }
                        }
                        stage("KEGG") {
                            when {
                                anyOf {
                                    expression { env.RUN_ALL != null }
                                    expression { env.KEGG != null }
                                }
                            }
                            steps {
                                sh '''
                                    SOURCE=kegg
                                    $DIPPER --sources $SOURCE
                                    scp ./out/${SOURCE}.ttl ./out/${SOURCE}_dataset.ttl $MONARCH_DATA_DEST
                                '''
                            }
                        }
                        stage("OMIA") {
                            when {
                                anyOf {
                                    expression { env.RUN_ALL != null }
                                    expression { env.OMIA != null }
                                }
                            }
                            steps {
                                sh '''
                                    SOURCE=omia
                                    $DIPPER --sources $SOURCE
                                    scp ./out/${SOURCE}.ttl ./out/${SOURCE}_dataset.ttl $MONARCH_DATA_DEST
                                '''
                            }
                        }
                        stage("Gene Reviews") {
                            when {
                                anyOf {
                                    expression { env.RUN_ALL != null }
                                    expression { env.GENE_REVIEWS != null }
                                }
                            }
                            steps {
                                dir('./data-boutique') {
                                    git(
                                        url: 'https://github.com/monarch-initiative/data-boutique.git',
                                        credentialsId: '3ca28d15-5fa8-46b1-a2ac-a5a483694f5b',
                                        branch: 'master'
                                    )
                                    sh '''
                                        SOURCE=genereviews
                                        cd .. && mkdir -p raw/genereviews/books
                                        cp ./data-boutique/GeneReviewsBooks/* ./raw/genereviews/books/
                                        $DIPPER --sources $SOURCE
                                    scp ./out/${SOURCE}.ttl ./out/${SOURCE}_dataset.ttl $MONARCH_DATA_DEST
                                    '''
                                }
                            }
                        }
                    }
                }
                stage("StringDb") {
                    when {
                        anyOf {
                            expression { env.RUN_ALL != null }
                            expression { env.STRING_DB != null }
                        }
                    }
                    steps {
                        sh '''
                            SOURCE=stringdb
                            $DIPPER --sources $SOURCE --taxon $COMMON_TAXON,10116 --version 11.0
                            scp ./out/string.ttl ./out/string_dataset.ttl $MONARCH_DATA_DEST
                        '''
                    }
                }
                stage("Panther") {
                    when {
                        anyOf {
                            expression { env.RUN_ALL != null }
                            expression { env.PANTHER != null }
                        }
                    }
                    steps {
                        sh '''
                            SOURCE=panther
                            mkdir -p raw/panther
                            $DIPPER --sources $SOURCE --taxon $COMMON_TAXON,10116,9913,9031,9796,9823,9615 --dest_fmt nt
                            scp ./out/${SOURCE}.nt ./out/${SOURCE}_dataset.ttl $MONARCH_DATA_DEST
                        '''
                    }
                }
                stage("AnimalQTLdb") {
                    when {
                        anyOf {
                            expression { env.RUN_ALL != null }
                            expression { env.ANIMALQTLDB != null }
                        }
                    }
                    steps {
                        sh '''
                            SOURCE=animalqtldb
                            $DIPPER --sources $SOURCE
                            scp ./out/${SOURCE}.ttl ./out/${SOURCE}_dataset.ttl $MONARCH_DATA_DEST
                        '''
                    }
                }
                stage("Bgee") {
                    when {
                        anyOf {
                            expression { env.RUN_ALL != null }
                            expression { env.BGEE != null }
                        }
                    }
                    steps {
                        sh '''
                            SOURCE=bgee
                            $DIPPER --sources $SOURCE --taxon $COMMON_TAXON,10116 --version bgee_v14_0 --limit 20

                            echo "check statement count and if well-formed?"
                            rapper -i turtle -c ./out/bgee.ttl
                            scp ./out/${SOURCE}.ttl ./out/${SOURCE}_dataset.ttl $MONARCH_DATA_DEST
                        '''
                    }
                }
                stage("FlyBase") {
                    when {
                        anyOf {
                            expression { env.RUN_ALL != null }
                            expression { env.FLYBASE != null }
                        }
                    }
                    steps {
                        sh '''
                            SOURCE=flybase
                            $DIPPER --sources $SOURCE
                            scp ./out/${SOURCE}.ttl ./out/${SOURCE}_dataset.ttl $MONARCH_DATA_DEST
                        '''
                    }
                }
                stage("Biogrid") {
                    when {
                        anyOf {
                            expression { env.RUN_ALL != null }
                            expression { env.BIOGRID != null }
                        }
                    }
                    steps {
                        sh '''
                            SOURCE=biogrid
                            $DIPPER --sources $SOURCE --taxon $COMMON_TAXON
                            scp ./out/${SOURCE}.ttl ./out/${SOURCE}_dataset.ttl $MONARCH_DATA_DEST
                        '''
                    }
                }
                stage("ClinVar") {
                    when {
                        anyOf {
                            expression { env.RUN_ALL != null }
                            expression { env.CLINVAR != null }
                        }
                    }
                    steps {
                        sh '''
                            mkdir -p out
                            mkdir -p raw/clinvar && cd raw/clinvar
                            # these are available via http in DipperCache too
                            # ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/ClinVarFullRelease_00-latest.xml.gz
                            wget --quiet --timestamping "$DIPPERCACHE/clinvar/ClinVarFullRelease_00-latest.xml.gz"
                            wget --quiet --timestamping "$DIPPERCACHE/clinvar/gene_condition_source_id"
                            cd -

                            export PYTHONPATH=.:$PYTHONPATH
                            venv/bin/python ./dipper/sources/ClinVar.py
                            scp ./out/clinvar.nt $MONARCH_DATA_DEST
                        '''
                    }
                }
                stage("Coriell") {
                    when {
                        anyOf {
                            expression { env.RUN_ALL != null }
                            expression { env.CORIELL != null }
                        }
                    }
                    steps {
                        sh '''
                            SOURCE=coriell
                            $DIPPER --sources $SOURCE
                            scp ./out/${SOURCE}.ttl ./out/${SOURCE}_dataset.ttl $MONARCH_DATA_DEST
                        '''
                    }
                }
                stage("CTD") {
                    when {
                        anyOf {
                            expression { env.RUN_ALL != null }
                            expression { env.CTD != null }
                        }
                    }
                    steps {
                        sh '''
                            SOURCE=ctd
                            $DIPPER --sources $SOURCE
                            scp ./out/${SOURCE}.ttl ./out/${SOURCE}_dataset.ttl $MONARCH_DATA_DEST
                        '''
                    }
                }
                stage("Ensembl") {
                    when {
                        anyOf {
                            expression { env.RUN_ALL != null }
                            expression { env.ENSEMBL != null }
                        }
                    }
                    steps {
                        sh '''
                            SOURCE=ensembl
                            $DIPPER --sources $SOURCE --taxon $COMMON_TAXON
                            scp ./out/${SOURCE}.ttl ./out/${SOURCE}_dataset.ttl $MONARCH_DATA_DEST
                        '''
                    }
                }
                stage("Elements of Morphology") {
                    when {
                        anyOf {
                            expression { env.RUN_ALL != null }
                            expression { env.EOM != null }
                        }
                    }
                    steps {
                        sh '''
                            SOURCE=eom
                            $DIPPER --sources $SOURCE
                            scp ./out/${SOURCE}.ttl ./out/${SOURCE}_dataset.ttl $MONARCH_DATA_DEST
                        '''
                    }
                }
                stage("Gene Ontology Associations") {
                    when {
                        anyOf {
                            expression { env.RUN_ALL != null }
                            expression { env.GO != null }
                        }
                    }
                    steps {
                        sh '''
                            SOURCE=go
                            $DIPPER --sources $SOURCE --dest_fmt nt --taxon \
                                $COMMON_TAXON,10116,4896,5052,559292,5782,9031,9615,9823,9913
                            scp ./out/${SOURCE}.nt ./out/${SOURCE}_dataset.ttl $MONARCH_DATA_DEST
                        '''
                    }
                }
                stage("GWAS Catalog") {
                    when {
                        anyOf {
                            expression { env.RUN_ALL != null }
                            expression { env.GWASCATALOG != null }
                        }
                    }
                    steps {
                        sh '''
                            SOURCE=gwascatalog
                            $DIPPER --sources $SOURCE
                            scp ./out/${SOURCE}.ttl ./out/${SOURCE}_dataset.ttl $MONARCH_DATA_DEST
                        '''
                    }
                }
                stage("HPO Annotations") {
                    when {
                        anyOf {
                            expression { env.RUN_ALL != null }
                            expression { env.HPOA != null }
                        }
                    }
                    steps {
                        sh '''
                            SOURCE=hpoa
                            $DIPPER --sources $SOURCE
                            scp ./out/${SOURCE}.ttl ./out/${SOURCE}_dataset.ttl $MONARCH_DATA_DEST
                        '''
                    }
                }
                stage("IMPC") {
                    when {
                        anyOf {
                            expression { env.RUN_ALL != null }
                            expression { env.IMPC != null }
                        }
                    }
                    steps {
                        sh '''
                            SOURCE=impc
                            $DIPPER --sources $SOURCE
                            scp ./out/${SOURCE}.ttl ./out/${SOURCE}_dataset.ttl $MONARCH_DATA_DEST
                        '''
                    }
                }
                stage("MGISlim") {
                    when {
                        anyOf {
                            expression { env.RUN_ALL != null }
                            expression { env.MGISLIM != null }
                        }
                    }
                    steps {
                        sh '''
                            SOURCE=mgislim
                            $DIPPER --sources $SOURCE
                            scp ./out/mgislim.ttl ./out/mgislim_dataset.ttl $MONARCH_DATA_DEST
                        '''
                    }
                }
                stage("MGI") {
                    when {
                        anyOf {
                            expression { env.RUN_ALL != null }
                            expression { env.MGI != null }
                        }
                    }
                    steps {
                        sh '''
                            SOURCE=mgi
                            $DIPPER --sources $SOURCE
                            scp ./out/${SOURCE}.ttl ./out/${SOURCE}_dataset.ttl $MONARCH_DATA_DEST
                        '''
                    }
                }
                stage("MMRRC") {
                    when {
                        anyOf {
                            expression { env.RUN_ALL != null }
                            expression { env.MMRRC != null }
                        }
                    }
                    steps {
                        sh '''
                            SOURCE=mmrrc
                            $DIPPER --sources $SOURCE
                            scp ./out/${SOURCE}.ttl ./out/${SOURCE}_dataset.ttl $MONARCH_DATA_DEST
                        '''
                    }
                }
                stage("Monarch Boutique") {
                    when {
                        anyOf {
                            expression { env.RUN_ALL != null }
                            expression { env.MONARCH != null }
                        }
                    }
                    steps {
                        dir('./data-boutique-b') {
                            git(
                                url: 'https://github.com/monarch-initiative/data-boutique.git',
                                credentialsId: '3ca28d15-5fa8-46b1-a2ac-a5a483694f5b',
                                branch: 'master'
                            )
                            sh '''
                                SOURCE=monarch
                                cd .. && mkdir -p raw/monarch
                                cp -r data-boutique-b/OMIA-disease-phenotype ./raw/monarch/
                                $DIPPER --sources $SOURCE
                                scp ./out/${SOURCE}.ttl ./out/${SOURCE}_dataset.ttl $MONARCH_DATA_DEST
                            '''
                        }
                    }
                }
                stage("Monochrom") {
                    when {
                        anyOf {
                            expression { env.RUN_ALL != null }
                            expression { env.MONOCHROM != null }
                        }
                    }
                    steps {
                        sh '''
                            SOURCE=monochrom
                            $DIPPER --sources $SOURCE --use_bnodes
                            scp ./out/${SOURCE}.ttl ./out/${SOURCE}_dataset.ttl $MONARCH_DATA_DEST
                        '''
                    }
                }
                stage("MPD") {
                    when {
                        anyOf {
                            expression { env.RUN_ALL != null }
                            expression { env.MPD != null }
                        }
                    }
                    steps {
                        sh '''
                            SOURCE=mpd
                            $DIPPER --sources $SOURCE
                            scp ./out/${SOURCE}.ttl ./out/${SOURCE}_dataset.ttl $MONARCH_DATA_DEST
                        '''
                    }
                }
                stage("Orphanet") {
                    when {
                        anyOf {
                            expression { env.RUN_ALL != null }
                            expression { env.ORPHANET != null }
                        }
                    }
                    steps {
                        sh '''
                            SOURCE=orphanet
                            $DIPPER --sources $SOURCE
                            scp ./out/${SOURCE}.ttl ./out/${SOURCE}_dataset.ttl $MONARCH_DATA_DEST
                        '''
                    }
                }
                stage("Reactome") {
                    when {
                        anyOf {
                            expression { env.RUN_ALL != null }
                            expression { env.REACTOME != null }
                        }
                    }
                    steps {
                        sh '''
                            SOURCE=reactome
                            $DIPPER --sources $SOURCE
                            scp ./out/${SOURCE}.ttl ./out/${SOURCE}_dataset.ttl $MONARCH_DATA_DEST
                        '''
                    }
                }
                stage("Wormbase") {
                    when {
                        anyOf {
                            expression { env.RUN_ALL != null }
                            expression { env.WORMBASE != null }
                        }
                    }
                    steps {
                        sh '''
                            SOURCE=wormbase
                            $DIPPER --sources $SOURCE --dest_fmt nt
                            scp ./out/${SOURCE}.nt ./out/${SOURCE}_dataset.ttl $MONARCH_DATA_DEST
                        '''
                    }
                }
                stage('Danio'){ /* the slim after the main ingest to use same files */
                    stages {
                        stage("ZFIN") {
                            when {
                                anyOf {
                                    expression { env.RUN_ALL != null }
                                    expression { env.ZFIN != null }
                                }
                            }
                            steps {
                                sh '''
                                    SOURCE=zfin
                                    $DIPPER --sources $SOURCE
                                    scp ./out/${SOURCE}.ttl ./out/${SOURCE}_dataset.ttl $MONARCH_DATA_DEST
                                '''
                            }
                        }
                        stage("ZFINSlim") {
                            when {
                                anyOf {
                                    expression { env.RUN_ALL != null }
                                    expression { env.ZFINSLIM != null }
                                }
                            }
                            steps {
                                sh '''
                                    SOURCE=zfinslim
                                    $DIPPER --sources $SOURCE
                                    scp ./out/${SOURCE}.ttl ./out/${SOURCE}_dataset.ttl $MONARCH_DATA_DEST
                                '''
                            }
                        }
                    }
                }
                stage("RGD") {
                    when {
                        anyOf {
                            expression { env.RUN_ALL != null }
                            expression { env.RGD != null }
                        }
                    }
                    steps {
                        sh '''
                            SOURCE=rgd
                            $DIPPER --sources $SOURCE
                            scp ./out/${SOURCE}.ttl ./out/${SOURCE}_dataset.ttl $MONARCH_DATA_DEST
                        '''
                    }
                }
                stage("SGD") {
                    when {
                        anyOf {
                            expression { env.RUN_ALL != null }
                            expression { env.SGD != null }
                        }
                    }
                    steps {
                        sh '''
                            SOURCE=sgd
                            $DIPPER --sources $SOURCE
                            scp ./out/${SOURCE}.ttl ./out/${SOURCE}_dataset.ttl $MONARCH_DATA_DEST
                        '''
                    }
                }
                /*
                stage("MyChem Info") {
                    when {
                        anyOf {
                            expression { env.RUN_ALL != null }
                            expression { env.MYCHEM != null }
                        }
                    }
                    steps {
                        sh '''
                            SOURCE=mychem
                            $DIPPER --sources $SOURCE
                            scp ./out/${SOURCE}.ttl ./out/${SOURCE}_dataset.ttl $MONARCH_DATA_DEST
                        '''
                    }
                }
                */
                stage("UCSCBands") {
                    when {
                        anyOf {
                            expression { env.RUN_ALL != null }
                            expression { env.UCSCBANDS != null }
                        }
                    }
                    steps {
                        sh '''
                            SOURCE=ucscbands
                            $DIPPER --sources $SOURCE
                            scp ./out/${SOURCE}.ttl ./out/${SOURCE}_dataset.ttl $MONARCH_DATA_DEST
                        '''
                    }
                }
                stage("Xenbase") {
                    when {
                        anyOf {
                            expression { env.RUN_ALL != null }
                            expression { env.XENBASE != null }
                        }
                    }
                    steps {
                        sh '''
                            SOURCE=xenbase
                            $DIPPER --sources $SOURCE
                            scp ./out/${SOURCE}.ttl ./out/${SOURCE}_dataset.ttl $MONARCH_DATA_DEST
                        '''
                    }
                }
            }
        }
        stage('Estatic'){
            when {
                expression { env.RELEASE != null }
            }
            steps {
                sh '''
                    echo "Taking a moment of silent contemplation ..."
                    sleep 1m
                    echo "relax perms"
                    chmod g+rw out/*
                    # Move Data to Monarch Archive
                    ./scripts/mdma.sh
                '''
            }
        }
        stage("Generate monarch merged owl file") {
            when {
                anyOf {
                    expression { env.RUN_ALL != null }
                    expression { env.MONARCH_OWL != null }
                }
            }
            steps {
                dir('./create-monarch-owl') {deleteDir()}
                dir('./create-monarch-owl') {
                    sh '''
                        wget --quiet --timestamping http://current.geneontology.org/bin/owltools
                        chmod +x owltools
                        java -Xmx100g -jar owltools http://purl.obolibrary.org/obo/upheno/monarch.owl --merge-import-closure --remove-disjoints --remove-equivalent-to-nothing-axioms -o monarch-merged.owl

                        # Hack to resolve https://github.com/monarch-initiative/monarch-ontology/issues/16
                        # Hack to normalize omim and hgnc IRIs

                        sed -i "/owl#ReflexiveProperty/d;\
                            s~http://purl.obolibrary.org/obo/OMIMPS_~http://www.omim.org/phenotypicSeries/PS~;\
                            s~http://purl.obolibrary.org/obo/OMIM_~http://omim.org/entry/~;\
                            s~http://identifiers.org/omim/~http://omim.org/entry/~;\
                            s~http://identifiers.org/hgnc/~https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:~;\
                            s~http://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=~https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:~;\
                            s~http://www.informatics.jax.org/marker/MGI:~http://www.informatics.jax.org/accession/MGI:~;\
                            s~http://www.ncbi.nlm.nih.gov/gene/~https://www.ncbi.nlm.nih.gov/gene/~; \
                            s~http://purl.obolibrary.org/obo/MESH_~http://id.nlm.nih.gov/mesh/~" \
                            ./monarch-merged.owl

                        scp monarch-merged.owl $MONARCH_DATA_FS:/var/www/data/owl/
                    '''
                }
            }
        }
    }
    /* alternative  to stage 'Estatic'
    post {
        // regression, aborted, failure, success, unstable, unsuccessful, cleanup //
        success {
            sh '''
                echo "Taking a moment of silent contemplation ..."
                sleep 1m
                echo "relax perms"
                chmod g+rw out/*
                # Move Data to Monarch Archive
                ./scripts/mdma.sh
            '''
        }
        cleanup {}
    }
    */
}
