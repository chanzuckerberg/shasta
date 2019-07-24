# Definitions
git_commit ?= $(shell git log --pretty=oneline -n 1 | cut -f1 -d " ")
name = tpesout/shasta
version = 1.0.0
tag = ${version}--${git_commit}

# Steps
build:
	# do the docker build
	docker build -t ${name}:${tag} .
	docker tag ${name}:${tag} ${name}:latest

push: build
	# Requires ~/.dockercfg
	docker push ${name}:${tag}
	docker push ${name}:latest
