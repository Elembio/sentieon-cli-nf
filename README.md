## sentieon-nf

Build+push docker

```
docker build -t sentieon-cli:202308.03 .sentieon-cli:202308.03
docker tag sentieon-cli:202308.03 147637153406.dkr.ecr.us-west-2.amazonaws.com/modules/sentieon-cli:202308.03
docker push 147637153406.dkr.ecr.us-west-2.amazonaws.com/modules/sentieon-cli:202308.03
```

Nextflow wrapper for Element Sentieon-NF (sentieon-cli)

```nextflow run . -profile docker,test```
