
container=hzsvela_maxflow

dockerfolder=dockerfiles

command="docker build $dockerfolder -t $container"

echo $command
$command

wd=$(pwd)

command="docker run --rm -it \
        --mount source=$wd/code,target=/code,type=bind \
        --mount source=$wd/input,target=/input,type=bind \
	$container /bin/bash"
#        --mount source=$wd/graphblast,target=/graphblast,type=bind \
#        --mount source=$wd/GraphBLAS,target=/GraphBLAS,type=bind \

#--mount source=$wd/gbtl,target=/gbtl,type=bind \

echo $command
$command
