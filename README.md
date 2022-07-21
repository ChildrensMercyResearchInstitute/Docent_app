#Docent app
##Build and push docker image
sh docker_build.sh 

##run the app
docker run -p 3838:3838 docker.io/man4ish/docent:latest

##browse the application
http://0.0.0.0:3838
