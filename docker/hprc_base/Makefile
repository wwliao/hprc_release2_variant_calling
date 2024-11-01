# Definitions
repository = wwliao
identifier = hprc_base
version = r2
git_commit ?= $(shell git log --pretty=format:"%h" -n 1)  # short commit hash
name = ${repository}/${identifier}
tag = ${version}-${git_commit}

# Steps
build:
	docker build -t ${name}:${tag} .
	docker tag ${name}:${tag} ${name}:latest

push:
	docker push ${name}:${tag}
	docker push ${name}:latest
