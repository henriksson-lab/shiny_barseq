docker:
	docker build -t shiny_barseq .
docker_testrun:
	docker run --rm -p 3838:3838 shiny_barseq
