# Definitions
git_commit = $(shell git log --pretty=oneline -n 1 | cut -f1 -d " ")
image_name = kishwar/shasta
shasta_version = 0.3.0
tag = ${shasta_version}--${git_commit}

# Steps
build:
	# do the docker build
	docker build -t ${image_name}:${tag} --build-arg git_commit=${git_commit} .
	docker tag ${image_name}:${tag} ${image_name}:latest

push: build
	# Requires ~/.dockercfg
	docker push ${image_name}:${tag}
	docker push ${image_name}:latest