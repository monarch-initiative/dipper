pipeline {

    agent any

    triggers {
        // Run every Monday at 5pm
        cron('H 17 * * 1')
    }

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

        MONARCH_DATA_FS = 'monarch-ttl-prod'
        DIPPER = 'venv/bin/python dipper-etl.py'
    }

    options {
        buildDiscarder(logRotator(numToKeepStr: '14'))
    }

    stages {

        stage('Build dipper package') {
            steps {
                git(
                    url: 'https://github.com/monarch-initiative/dipper.git',
                    branch: 'master'
                )
                dir('./config') {
                    git(
                        url: 'https://github.com/monarch-initiative/configs.git',
                        credentialsId: '3ca28d15-5fa8-46b1-a2ac-a5a483694f5b',
                        branch: 'master'
                    )
                    sh '''
                        cd .. && cp config/Dipper/conf.yaml ./dipper/conf.yaml
                        virtualenv -p /usr/bin/python3 venv
                        venv/bin/pip install -r requirements.txt
                        venv/bin/pip install -r requirements/all-sources.txt
                        
                        # Clean up previous runs
                        sudo rm -rf ./out/
                        sudo rm -rf ./raw/
                    '''
                }
            }
        }
        stage('Generate monarch owl and rdf') {
            steps {
                parallel(
                    "Generate monarch merged owl file": {
                        dir('./create-monarch-owl') {
                            sh """
                                wget http://build.berkeleybop.org/job/owltools/lastSuccessfulBuild/artifact/OWLTools-Runner/target/owltools

                                chmod +x owltools

                                java -Xmx100g -jar owltools http://purl.obolibrary.org/obo/upheno/monarch.owl --merge-import-closure --remove-disjoints --remove-equivalent-to-nothing-axioms -o monarch-merged.owl

                                # Hack to resolve https://github.com/monarch-initiative/monarch-ontology/issues/16
                                sed -i "/owl#ReflexiveProperty/d" ./monarch-merged.owl

                                # Hack to normalize omim and hgnc IRIs
                                # https://github.com/monarch-initiative/dipper/issues/700
                                sed -i 's/http:\\/\\/purl.obolibrary.org\\/obo\\/OMIMPS_/http:\\/\\/www.omim.org\\/phenotypicSeries\\//' ./monarch-merged.owl
                                sed -i 's/http:\\/\\/purl.obolibrary.org\\/obo\\/OMIM_/http:\\/\\/omim.org\\/entry\\//' ./monarch-merged.owl
                                sed -i 's/http:\\/\\/identifiers.org\\/omim\\//http:\\/\\/omim.org\\/entry\\//' ./monarch-merged.owl
                                sed -i 's/http:\\/\\/identifiers.org\\/hgnc\\//http:\\/\\/www.genenames.org\\/cgi-bin\\/gene_symbol_report?hgnc_id=/' ./monarch-merged.owl

                                scp monarch-merged.owl monarch@$MONARCH_DATA_FS:/var/www/data/owl/
                            """
                        }
                    },
                    "ETL StringDb": {
                        sh '''
                            $DIPPER --sources stringdb --taxon 6239,9606,10090,7955,7227,10116 --version 11.0
                            head -1000 ./out/string.ttl > ./out/string_test.ttl
                        '''
                    },
                    "ETL Panther": {
                        sh '''
                            $DIPPER --sources panther --taxon 9913,6239,9031,7955,7227,9796,9606,10090,9823,10116,8364,9615
                        '''
                    },
                    "ETL AnimalQTLdb": {
                        sh '''
                            $DIPPER --sources animalqtldb
                        '''
                    },
                    "ETL Bgee": {
                        sh '''
                            # Test run
                            $DIPPER --sources bgee --limit 1 --taxon 9606 --version bgee_v13_2
                            mv ./out/bgee.ttl ./out/bgee_test.ttl

                            $DIPPER --sources bgee --limit 20 --taxon 9606,10090,7227,6239,7955,10116 # --version bgee_v13_2

                            echo "check statement count and if well-formed?"
                            rapper -i turtle -c ./out/bgee.ttl
                        '''
                    },
                    "ETL FlyBase": {
                        sh '''
                            $DIPPER --sources flybase --dest_fmt nt
                        '''
                    },
                    "ETL Biogrid": {
                        sh '''
                            $DIPPER --sources biogrid --taxon 9606,10090,7955,7227,6239
                        '''
                    },
                    "ETL ClinVar": {
                        sh '''
                            mkdir -p out
                            mkdir -p raw && cd raw
                            mkdir -p clinvarxml_alpha && cd clinvarxml_alpha
                            wget -q --timestamping ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/ClinVarFullRelease_00-latest.xml.gz
                            wget -q --timestamping ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/gene_condition_source_id
                            cd ../..

                            export PYTHONPATH=.:$PYTHONPATH
                            venv/bin/python ./dipper/sources/ClinVarXML_alpha.py
                        '''
                    },
                    "ETL Coriell": {
                        sh '''
                            $DIPPER --sources coriell
                        '''
                    },
                    "ETL CTD": {
                        sh '''
                            $DIPPER --sources ctd
                        '''
                    },
                    "ETL Ensembl": {
                        sh '''
                            $DIPPER --sources ensembl --taxon 9606,10090,7955,7227,6239
                        '''
                    },
                    "ETL Elements of Morphology": {
                        sh '''
                            $DIPPER --sources eom
                        '''
                    },
                    "ETL Gene Reviews": {
                        dir('./data-boutique') {
                            git(
                                url: 'https://github.com/monarch-initiative/data-boutique.git',
                                credentialsId: '3ca28d15-5fa8-46b1-a2ac-a5a483694f5b',
                                branch: 'master'
                            )
                            sh '''
                                cd .. && mkdir -p raw/genereviews/books
                                cp ./data-boutique/GeneReviewsBooks/* ./raw/genereviews/books/
                                $DIPPER --sources genereviews
                            '''
                        }
                    },
                    "ETL Gene Ontology Associations": {
                        sh '''
                            $DIPPER --sources go --taxon \
                                10090,10116,4896,5052,559292,5782,6239,7227,7955,9031,9606,9615,9823,9913
                        '''
                    },
                    "ETL GWAS Catalog": {
                        sh '''
                            $DIPPER --sources gwascatalog --skip_tests   
                        '''
                    },
                    "ETL HGNC": {
                        sh '''
                            $DIPPER --sources hgnc
                        '''
                    },
                    "ETL HPO Annotations": {
                        sh '''
                            $DIPPER --sources hpoa  
                        '''
                    },
                    "ETL IMPC": {
                        sh '''
                            $DIPPER --sources impc
                        '''
                    },
                    "ETL KEGG": {
                        sh '''
                            $DIPPER --sources kegg      
                        '''
                    },
                    "ETL MGI Slim": {
                        sh '''
                            $DIPPER --sources mgi-slim
                        '''
                    },
                    "ETL MGI": {
                        sh '''
                            $DIPPER --sources mgi --skip_tests
                        '''
                    },
                    "ETL MMRRC": {
                        sh '''
                            $DIPPER --sources mmrrc      
                        '''
                    },
                    "ETL Monarch Boutique": {
                        dir('./data-boutique-b') {
                            git(
                                url: 'https://github.com/monarch-initiative/data-boutique.git',
                                credentialsId: '3ca28d15-5fa8-46b1-a2ac-a5a483694f5b',
                                branch: 'master'
                            )
                            sh '''
                                cd .. && mkdir -p raw/monarch
                                cp -r data-boutique-b/OMIA-disease-phenotype ./raw/monarch/
                                $DIPPER --sources monarch --skip_tests
                            '''
                        }
                    },
                    "ETL monochrom": {
                        sh '''
                            $DIPPER --sources monochrom          
                        '''
                    },
                    "ETL MPD": {
                        sh '''
                            $DIPPER --sources mpd --skip_tests    
                        '''
                    },
                    "ETL NCBI Gene": {
                        sh ''' 
                            $DIPPER --sources ncbigene --taxon \
                            28377,3702,9913,6239,9615,9031,7955,44689,7227,9796,9606,9544,13616,10090,9258,9598,9823,10116,4896,31033,8364,9685,559292
                        '''
                    },
                    "ETL OMIM": {
                        sh '''
                            $DIPPER --sources omim -q
                        '''
                    },
                    "ETL Orphanet": {
                        sh '''    
                            $DIPPER --sources orphanet        
                        '''
                    },
                    "ETL Reactome": {
                        sh '''
                            $DIPPER --sources reactome     
                        '''
                    },
                    "ETL Wormbase": {
                        sh '''
                            $DIPPER --sources wormbase --dest_fmt nt       
                        '''
                    },
                    "ETL ZFIN Slim": {
                        sh '''
                            $DIPPER --sources zfinslim
                        '''
                    },
                    "ETL ZFIN": {
                        sh '''
                            $DIPPER --sources zfin
                        '''
                    },
                    "ETL RGD": {
                        sh '''
                            $DIPPER --sources rgd   
                        '''
                    },
                    "ETL SGD": {
                        sh '''
                            $DIPPER --sources sgd
                        '''
                    },
                    "ETL MyChem Info": {
                        sh '''
                            $DIPPER --sources mychem
                        '''
                    },
                    "ETL OMIA": {
                        sh '''
                            $DIPPER --sources omia
                        '''
                    },
                    "ETL UCSCBands": {
                        sh '''
                            $DIPPER --sources ucscbands
                        '''
                    },
                )
            }
        }
    }

    post {
        always {
            // Copy files
            sh '''
                scp ./out/* monarch@$MONARCH_DATA_FS:/var/www/data/ttl/
            '''
        }
    }
}
