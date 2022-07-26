#Docent app
##Build and push docker image
sh docker_build.sh

##run the app
docker run -p 3838:3838 -v /path to data_folder/pbmc_2018-05-24_dev:/pbmc_2018-05-24_dev docker.io/man4ish/docent:latest

### path to data folder on gpfs
/home/mkumar1/app_backup/backup_20220708/Docent_app/pbmc_2018-05-24_dev

### Example
docker run -p 3838:3838 -v /home/mkumar1/app_backup/backup_20220708/Docent_app/pbmc_2018-05-24_dev:/pbmc_2018-05-24_dev docker.io/man4ish/docent:latest


##browse the application
http://0.0.0.0:3838
