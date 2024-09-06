check_build_push_image () (
  set -euo pipefail
  IMAGE="$1"
  CONTEXT="$2"
  if ! DOCKER_CLI_EXPERIMENTAL=enabled docker manifest inspect ${IMAGE} >/dev/null; then docker build -t ${IMAGE} ${CONTEXT} && docker push ${IMAGE}; fi
)

push_semver_tags () (
  set -euo pipefail
  docker push ${CI_REGISTRY_IMAGE}:${CI_COMMIT_TAG}
  major=`echo ${CI_COMMIT_TAG} | cut -d . -f 1`
  minor=`echo ${CI_COMMIT_TAG} | cut -d . -f 2`
  patch=`echo ${CI_COMMIT_TAG} | cut -d . -f 3`
  docker tag ${CI_REGISTRY_IMAGE}:${major}.${minor}.${patch} ${CI_REGISTRY_IMAGE}:${major}.${minor}
  docker tag ${CI_REGISTRY_IMAGE}:${major}.${minor}.${patch} ${CI_REGISTRY_IMAGE}:${major} 
  docker push ${CI_REGISTRY_IMAGE}:${major}.${minor}
  docker push ${CI_REGISTRY_IMAGE}:${major}
)
