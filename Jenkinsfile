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
    }

    options {
        buildDiscarder(logRotator(numToKeepStr: '14'))
    }

    stages {

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

                                # Hack to normalize omim IRIs
                                # https://github.com/monarch-initiative/dipper/issues/700
                                sed -i 's/http:\\/\\/purl.obolibrary.org\\/obo\\/OMIMPS_/http:\\/\\/www.omim.org\\/phenotypicSeries\\//' ./monarch-merged.owl
                                sed -i 's/http:\\/\\/purl.obolibrary.org\\/obo\\/OMIM_/http:\\/\\/omim.org\\/entry\\//' ./monarch-merged.owl
                                sed -i 's/http:\\/\\/identifiers.org\\/omim\\//http:\\/\\/omim.org\\/entry\\//' ./monarch-merged.owl

                                scp monarch-merged.owl monarch@$MONARCH_DATA_FS:/var/www/data/owl/
                            """
                        }
                    },
                    "ETL StringDb": {
                        dir('./build-stringdb-rdf') {
                            git(
                                url: 'https://github.com/monarch-initiative/dipper.git',
                                branch: 'master'
                            )
                            sh '''
                                virtualenv -p /usr/bin/python3 venv
                                venv/bin/pip install -r requirements.txt
                                venv/bin/pip install -r requirements/all-sources.txt
                                venv/bin/python dipper-etl.py --sources stringdb --taxon 6239,9606,10090,7955,7227,10116 --version 11.0
                                head -1000 ./out/string.ttl > ./out/string_test.ttl
                                scp ./out/* monarch@$MONARCH_DATA_FS:/var/www/data/ttl/
                            '''
                        }
                    },
                    "ETL Panther": {
                        dir('./build-panther-rdf') {
                            git(
                                url: 'https://github.com/monarch-initiative/dipper.git',
                                branch: 'master'
                            )
                            sh '''
                                virtualenv -p /usr/bin/python3 venv
                                venv/bin/pip install -r requirements.txt
                                venv/bin/pip install -r requirements/all-sources.txt
                                venv/bin/python dipper-etl.py --sources panther --taxon 9913,6239,9031,7955,7227,9796,9606,10090,9823,10116,8364,9615
                                scp ./out/* monarch@$MONARCH_DATA_FS:/var/www/data/ttl/
                            '''
                        }
                    },
                    "ETL AnimalQTLdb": {
                        dir('./build-animalqtldb-rdf') {
                            git(
                                url: 'https://github.com/monarch-initiative/dipper.git',
                                branch: 'master'
                            )
                            sh '''
                                virtualenv -p /usr/bin/python3 venv
                                venv/bin/pip install -r requirements.txt
                                venv/bin/pip install -r requirements/all-sources.txt
                                venv/bin/python dipper-etl.py --sources animalqtldb
                                scp ./out/* monarch@$MONARCH_DATA_FS:/var/www/data/ttl/
                            '''
                        }
                    },
                    "ETL Bgee": {
                        dir('./build-bgee-rdf') {
                            git(
                                url: 'https://github.com/monarch-initiative/dipper.git',
                                branch: 'master'
                            )
                            sh '''
                                virtualenv -p /usr/bin/python3 venv
                                venv/bin/pip install -r requirements.txt
                                venv/bin/pip install -r requirements/all-sources.txt

                                # Test run
                                venv/bin/python dipper-etl.py --sources bgee --limit 1 --taxon 9606 --version bgee_v13_2
                                mv ./out/bgee.ttl ./out/bgee_test.ttl

                                venv/bin/python dipper-etl.py --sources bgee --limit 20 --taxon 9606,10090,7227,6239,7955,10116 # --version bgee_v13_2

                                echo "check statement count and if well-formed?"
                                rapper -i turtle -c ./out/bgee.ttl

                                scp ./out/* monarch@$MONARCH_DATA_FS:/var/www/data/ttl/
                            '''
                        }
                    },
                    "ETL FlyBase": {
                        dir('./build-flybase-rdf') {
                            git(
                                url: 'https://github.com/monarch-initiative/dipper.git',
                                branch: '41a5e2d'
                            )
                            sh '''
                                virtualenv -p /usr/bin/python3 venv
                                venv/bin/pip install -r requirements.txt
                                venv/bin/pip install -r requirements/all-sources.txt
                                venv/bin/python dipper-etl.py --sources flybase --dest_fmt nt
                                scp ./out/* monarch@$MONARCH_DATA_FS:/var/www/data/ttl/
                            '''
                        }
                    },
                    "ETL Biogrid": {
                        dir('./build-biogrid-rdf') {
                            git(
                                url: 'https://github.com/monarch-initiative/dipper.git',
                                branch: 'master'
                            )
                            sh '''
                                virtualenv -p /usr/bin/python3 venv
                                venv/bin/pip install -r requirements.txt
                                venv/bin/pip install -r requirements/all-sources.txt
                                venv/bin/python dipper-etl.py --sources biogrid --taxon 9606,10090,7955,7227,6239
                                scp ./out/* monarch@$MONARCH_DATA_FS:/var/www/data/ttl/
                            '''
                        }
                    },
                    "ETL ClinVar": {
                        dir('./build-clinvar-rdf') {
                            git(
                                url: 'https://github.com/monarch-initiative/dipper.git',
                                branch: 'master'
                            )
                            sh '''
                                virtualenv -p /usr/bin/python3 venv
                                venv/bin/pip install -r requirements.txt
                                venv/bin/pip install -r requirements/all-sources.txt

                                mkdir -p out
                                mkdir -p raw && cd raw
                                mkdir -p clinvarxml_alpha && cd clinvarxml_alpha
                                wget -q --timestamping ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/ClinVarFullRelease_00-latest.xml.gz
                                wget -q --timestamping ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/gene_condition_source_id
                                cd ../..

                                venv/bin/python ./dipper/sources/ClinVarXML_alpha.py
                                venv/bin/python ./scripts/add-properties2turtle.py --input ./out/clinvarxml_alpha.nt --output ./out/clinvar.ttl --input_format nt --output_format turtle

                                scp ./out/* monarch@$MONARCH_DATA_FS:/var/www/data/ttl/
                            '''
                        }
                    },
                    "ETL Coriell": {
                        dir('./build-coriell-rdf/') {
                            git(
                                url: 'https://github.com/monarch-initiative/dipper.git',
                                branch: 'master'
                            )
                        }
                        dir('./build-coriell-rdf/config') {
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
                                venv/bin/python dipper-etl.py --sources coriell

                                scp ./out/* monarch@$MONARCH_DATA_FS:/var/www/data/ttl/
                            '''
                        }
                    },
                    "ETL CTD": {
                        dir('./build-ctd-rdf') {
                            git(
                                url: 'https://github.com/monarch-initiative/dipper.git',
                                branch: 'master'
                            )
                            sh '''
                                virtualenv -p /usr/bin/python3 venv
                                venv/bin/pip install -r requirements.txt
                                venv/bin/pip install -r requirements/all-sources.txt
                                venv/bin/python dipper-etl.py --sources ctd
                                scp ./out/* monarch@$MONARCH_DATA_FS:/var/www/data/ttl/
                            '''
                        }
                    },
                    "ETL Ensembl": {
                        dir('./build-ensembl-rdf') {
                            git(
                                url: 'https://github.com/monarch-initiative/dipper.git',
                                branch: 'master'
                            )
                            sh '''
                                virtualenv -p /usr/bin/python3 venv
                                venv/bin/pip install -r requirements.txt
                                venv/bin/pip install -r requirements/all-sources.txt
                                venv/bin/python dipper-etl.py --sources ensembl --taxon 9606,10090,7955,7227,6239
                                scp ./out/* monarch@$MONARCH_DATA_FS:/var/www/data/ttl/
                            '''
                        }
                    },
                    "ETL Elements of Morphology": {
                        dir('./build-eom-rdf/') {
                            git(
                                url: 'https://github.com/monarch-initiative/dipper.git',
                                branch: 'master'
                            )
                        }
                        dir('./build-eom-rdf/config') {
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
                                venv/bin/python dipper-etl.py --sources eom

                                scp ./out/* monarch@$MONARCH_DATA_FS:/var/www/data/ttl/
                            '''
                        }
                    },
                    "ETL Gene Reviews": {
                        dir('./build-genereviews-rdf/') {
                            git(
                                url: 'https://github.com/monarch-initiative/dipper.git',
                                branch: 'master'
                            )
                        }
                        dir('./build-genereviews-rdf/config') {
                            git(
                                url: 'https://github.com/monarch-initiative/configs.git',
                                credentialsId: '3ca28d15-5fa8-46b1-a2ac-a5a483694f5b',
                                branch: 'master'
                            )
                            sh '''
                                cd .. && cp config/Dipper/conf.yaml ./dipper/conf.yaml
                            '''
                        }

                        dir('./build-genereviews-rdf/data-boutique') {
                            git(
                                url: 'https://github.com/monarch-initiative/data-boutique.git',
                                credentialsId: '3ca28d15-5fa8-46b1-a2ac-a5a483694f5b',
                                branch: 'master'
                            )
                            sh '''
                                cd .. && mkdir -p raw/genereviews/books
                                cp ./data-boutique/GeneReviewsBooks/* ./raw/genereviews/books/

                                virtualenv -p /usr/bin/python3 venv
                                venv/bin/pip install -r requirements.txt
                                venv/bin/pip install -r requirements/all-sources.txt
                                venv/bin/python dipper-etl.py --sources genereviews

                                scp ./out/* monarch@$MONARCH_DATA_FS:/var/www/data/ttl/
                            '''
                        }
                    },
                    "ETL Gene Ontology Associations": {
                        dir('./build-goa-rdf') {
                            git(
                                url: 'https://github.com/monarch-initiative/dipper.git',
                                branch: 'master'
                            )
                            sh '''
                                virtualenv -p /usr/bin/python3 venv
                                venv/bin/pip install -r requirements.txt
                                venv/bin/pip install -r requirements/all-sources.txt
                                venv/bin/python dipper-etl.py --sources go --taxon 7955,9606,10090,6239,7227,10116,559292,4896,4932
                                scp ./out/* monarch@$MONARCH_DATA_FS:/var/www/data/ttl/
                            '''
                        }
                    },
                    "ETL GWAS Catalog": {
                        dir('./build-gwascatalog-rdf') {
                            git(
                                url: 'https://github.com/monarch-initiative/dipper.git',
                                branch: 'master'
                            )
                            sh '''
                                virtualenv -p /usr/bin/python3 venv
                                venv/bin/pip install -r requirements.txt
                                venv/bin/pip install -r requirements/all-sources.txt
                                venv/bin/python dipper-etl.py --sources gwascatalog --skip_tests
                                scp ./out/* monarch@$MONARCH_DATA_FS:/var/www/data/ttl/
                            '''
                        }
                    },
                    "ETL HGNC": {
                        dir('./build-hgnc-rdf') {
                            git(
                                url: 'https://github.com/monarch-initiative/dipper.git',
                                branch: 'master'
                            )
                        }
                        dir('./build-hgnc-rdf/config') {
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
                                venv/bin/python dipper-etl.py --sources hgnc
                                scp ./out/* monarch@$MONARCH_DATA_FS:/var/www/data/ttl/
                            '''
                        }
                    },
                    "ETL HPO Annotations": {
                        dir('./build-hpoa-rdf/') {
                            git(
                                url: 'https://github.com/monarch-initiative/dipper.git',
                                branch: 'master'
                            )
                        }
                        dir('./build-hpoa-rdf/config') {
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
                                venv/bin/python dipper-etl.py --sources hpoa

                                scp ./out/* monarch@$MONARCH_DATA_FS:/var/www/data/ttl/
                            '''
                        }
                    },
                    "ETL IMPC": {
                        dir('./build-impc-rdf') {
                            git(
                                url: 'https://github.com/monarch-initiative/dipper.git',
                                branch: 'master'
                            )
                            sh '''
                                virtualenv -p /usr/bin/python3 venv
                                venv/bin/pip install -r requirements.txt
                                venv/bin/pip install -r requirements/all-sources.txt
                                venv/bin/python dipper-etl.py --sources impc
                                scp ./out/* monarch@$MONARCH_DATA_FS:/var/www/data/ttl/
                            '''
                        }
                    },
                    "ETL KEGG": {
                        dir('./build-kegg-rdf') {
                            git(
                                url: 'https://github.com/monarch-initiative/dipper.git',
                                branch: 'master'
                            )
                        }
                        dir('./build-krgg-rdf/config') {
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
                                venv/bin/python dipper-etl.py --sources kegg
                                scp ./out/* monarch@$MONARCH_DATA_FS:/var/www/data/ttl/
                            '''
                        }
                    },
                    "ETL MGI Slim": {
                        dir('./build-mgislim-rdf') {
                            git(
                                url: 'https://github.com/monarch-initiative/dipper.git',
                                branch: 'master'
                            )
                            sh '''
                                virtualenv -p /usr/bin/python3 venv
                                venv/bin/pip install -r requirements.txt
                                venv/bin/pip install -r requirements/all-sources.txt
                                venv/bin/python dipper-etl.py --sources mgi-slim
                                scp ./out/* monarch@$MONARCH_DATA_FS:/var/www/data/ttl/
                            '''
                        }
                    },
                    "ETL MGI": {
                        dir('./build-mgi-rdf/') {
                            git(
                                url: 'https://github.com/monarch-initiative/dipper.git',
                                branch: 'master'
                            )
                        }
                        dir('./build-mgi-rdf/config') {
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
                                venv/bin/python dipper-etl.py --sources mgi --skip_tests

                                scp ./out/* monarch@$MONARCH_DATA_FS:/var/www/data/ttl/
                            '''
                        }
                    },
                    "ETL MMRRC": {
                        dir('./build-mmrrc-rdf') {
                            git(
                                url: 'https://github.com/monarch-initiative/dipper.git',
                                branch: 'master'
                            )
                            sh '''
                                virtualenv -p /usr/bin/python3 venv
                                venv/bin/pip install -r requirements.txt
                                venv/bin/pip install -r requirements/all-sources.txt
                                venv/bin/python dipper-etl.py --sources mmrrc
                                scp ./out/* monarch@$MONARCH_DATA_FS:/var/www/data/ttl/
                            '''
                        }
                    },
                    "ETL Monarch Boutique": {
                        dir('./build-monarch-rdf/') {
                            git(
                                url: 'https://github.com/monarch-initiative/dipper.git',
                                branch: 'master'
                            )
                        }
                        dir('./build-monarch-rdf/data-boutique') {
                            git(
                                url: 'https://github.com/monarch-initiative/data-boutique.git',
                                credentialsId: '3ca28d15-5fa8-46b1-a2ac-a5a483694f5b',
                                branch: 'master'
                            )
                            sh '''
                                cd .. && mkdir -p raw/monarch
                                cp -r data-boutique/OMIA-disease-phenotype ./raw/monarch/

                                virtualenv -p /usr/bin/python3 venv
                                venv/bin/pip install -r requirements.txt
                                venv/bin/pip install -r requirements/all-sources.txt
                                venv/bin/python dipper-etl.py --sources monarch --skip_tests

                                scp ./out/* monarch@$MONARCH_DATA_FS:/var/www/data/ttl/
                            '''
                        }
                    },
                    "ETL monochrom": {
                        dir('./build-monochrom-rdf') {
                            git(
                                url: 'https://github.com/monarch-initiative/dipper.git',
                                branch: 'master'
                            )
                            sh '''
                                virtualenv -p /usr/bin/python3 venv
                                venv/bin/pip install -r requirements.txt
                                venv/bin/pip install -r requirements/all-sources.txt
                                venv/bin/python dipper-etl.py --sources monochrom
                                scp ./out/* monarch@$MONARCH_DATA_FS:/var/www/data/ttl/
                            '''
                        }
                    },
                    "ETL MPD": {
                        dir('./build-mpd-rdf') {
                            git(
                                url: 'https://github.com/monarch-initiative/dipper.git',
                                branch: 'master'
                            )
                            sh '''
                                virtualenv -p /usr/bin/python3 venv
                                venv/bin/pip install -r requirements.txt
                                venv/bin/pip install -r requirements/all-sources.txt
                                venv/bin/python dipper-etl.py --sources mpd --skip_tests
                                scp ./out/* monarch@$MONARCH_DATA_FS:/var/www/data/ttl/
                            '''
                        }
                    },
                    "ETL NCBI Gene": {
                        dir('./build-ncbigene-rdf') {
                            git(
                                url: 'https://github.com/monarch-initiative/dipper.git',
                                branch: 'master'
                            )
                        }
                        dir('./build-ncbigene-rdf/config') {
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
                                venv/bin/python dipper-etl.py --sources ncbigene --taxon 28377,3702,9913,6239,9615,9031,7955,44689,7227,9796,9606,9544,13616,10090,9258,9598,9823,10116,4896,31033,8364,4932,9685,559292
                                scp ./out/* monarch@$MONARCH_DATA_FS:/var/www/data/ttl/
                            '''
                        }
                    },
                    "ETL OMIM": {
                        dir('./build-omim-rdf/') {
                            git(
                                url: 'https://github.com/monarch-initiative/dipper.git',
                                branch: 'master'
                            )
                        }
                        dir('./build-omim-rdf/config') {
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
                                venv/bin/python dipper-etl.py --sources omim -q

                                scp ./out/omim.ttl monarch@$MONARCH_DATA_FS:/var/www/data/ttl/
                                scp ./out/omim_dataset.ttl monarch@$MONARCH_DATA_FS:/var/www/data/ttl/
                            '''
                        }
                    },
                    "ETL Orphanet": {
                        dir('./build-orphanet-rdf/') {
                            git(
                                url: 'https://github.com/monarch-initiative/dipper.git',
                                branch: 'master'
                            )
                            sh '''
                                virtualenv -p /usr/bin/python3 venv
                                venv/bin/pip install -r requirements.txt
                                venv/bin/pip install -r requirements/all-sources.txt
                                venv/bin/python dipper-etl.py --sources orphanet

                                scp ./out/* monarch@$MONARCH_DATA_FS:/var/www/data/ttl/
                            '''
                        }
                    },
                    "ETL Reactome": {
                        dir('./build-reactome-rdf/') {
                            git(
                                url: 'https://github.com/monarch-initiative/dipper.git',
                                branch: 'master'
                            )
                            sh '''
                                virtualenv -p /usr/bin/python3 venv
                                venv/bin/pip install -r requirements.txt
                                venv/bin/pip install -r requirements/all-sources.txt
                                venv/bin/python dipper-etl.py --sources reactome

                                scp ./out/* monarch@$MONARCH_DATA_FS:/var/www/data/ttl/
                            '''
                        }
                    },
                    "ETL Wormbase": {
                        dir('./build-wormbase-rdf/') {
                            git(
                                url: 'https://github.com/monarch-initiative/dipper.git',
                                branch: 'master'
                            )
                            sh '''
                                virtualenv -p /usr/bin/python3 venv
                                venv/bin/pip install -r requirements.txt
                                venv/bin/pip install -r requirements/all-sources.txt
                                venv/bin/python dipper-etl.py --sources wormbase --dest_fmt nt

                                scp ./out/* monarch@$MONARCH_DATA_FS:/var/www/data/ttl/
                            '''
                        }
                    },
                    "ETL ZFIN Slim": {
                        dir('./build-zfinslim-rdf/') {
                            git(
                                url: 'https://github.com/monarch-initiative/dipper.git',
                                branch: 'master'
                            )
                            sh '''
                                virtualenv -p /usr/bin/python3 venv
                                venv/bin/pip install -r requirements.txt
                                venv/bin/pip install -r requirements/all-sources.txt
                                venv/bin/python dipper-etl.py --sources zfinslim

                                scp ./out/* monarch@$MONARCH_DATA_FS:/var/www/data/ttl/
                            '''
                        }
                    },
                    "ETL ZFIN": {
                        dir('./build-zfin-rdf/') {
                            git(
                                url: 'https://github.com/monarch-initiative/dipper.git',
                                branch: 'master'
                            )
                            sh '''
                                virtualenv -p /usr/bin/python3 venv
                                venv/bin/pip install -r requirements.txt
                                venv/bin/pip install -r requirements/all-sources.txt
                                venv/bin/python dipper-etl.py --sources zfin

                                scp ./out/* monarch@$MONARCH_DATA_FS:/var/www/data/ttl/
                            '''
                        }
                    },
                    "ETL RGD": {
                        dir('./build-rgd-rdf/') {
                            git(
                                url: 'https://github.com/monarch-initiative/dipper.git',
                                branch: 'master'
                            )
                            sh '''
                                virtualenv -p /usr/bin/python3 venv
                                venv/bin/pip install -r requirements.txt
                                venv/bin/pip install -r requirements/all-sources.txt
                                venv/bin/python dipper-etl.py --sources rgd

                                scp ./out/* monarch@$MONARCH_DATA_FS:/var/www/data/ttl/
                            '''
                        }
                    },
                    "ETL SGD": {
                        dir('./build-sgd-rdf') {
                            git(
                                url: 'https://github.com/monarch-initiative/dipper.git',
                                branch: 'master'
                            )
                            sh '''
                                virtualenv -p /usr/bin/python3 venv
                                venv/bin/pip install -r requirements.txt
                                venv/bin/pip install -r requirements/all-sources.txt
                                venv/bin/python dipper-etl.py --sources sgd
                                scp ./out/* monarch@$MONARCH_DATA_FS:/var/www/data/ttl/
                            '''
                        }
                    },
                    "ETL MyChem Info": {
                        dir('./build-mychem-rdf/') {
                            git(
                                url: 'https://github.com/monarch-initiative/dipper.git',
                                branch: 'master'
                            )
                            sh '''
                                virtualenv -p /usr/bin/python3 venv
                                venv/bin/pip install -r requirements.txt
                                venv/bin/pip install -r requirements/all-sources.txt
                                venv/bin/python dipper-etl.py --sources mychem

                                scp ./out/* monarch@$MONARCH_DATA_FS:/var/www/data/ttl/
                            '''
                        }
                    },
                )
            }
        }
    }
    /*
    post {
        success {
        }
        changed {
        }
        failure {
        }
    }*/
}
