FROM node:latest as prepare_package_json
WORKDIR /app
COPY . /app
RUN find /app -type f -a \! \( -name "package.json" -o -name "package-lock.json" \) | xargs rm -f

FROM node:latest as prepare_requirements_txt
WORKDIR /app
COPY . /app
RUN find /app -type f -a \! \( -name "requirements.txt" -o -name "enumerate-requirements.ts" \) | xargs rm -f

FROM node:latest as prepare_npm_i
WORKDIR /app
COPY --from=prepare_package_json /app /app
RUN npm i

FROM node:latest as prepare_requirements_txt_complete
WORKDIR /app
COPY --from=prepare_requirements_txt /app /app
COPY --from=prepare_npm_i /app /app
RUN npm run codegen:requirements
RUN mv /app/requirements.txt /tmp/requirements.txt && rm -r /app && mkdir /app && mv /tmp/requirements.txt /app

FROM node:latest as prepare_build
WORKDIR /app
COPY --from=prepare_npm_i /app /app
COPY . /app
ENV LANDING_PAGE=/graph/extend
RUN npm run build

FROM node:latest as app
WORKDIR /app
ENV PYTHON_BIN="python3"
RUN set -x \
  && apt-get -y update \
  && apt-get -y install \
    python3-dev \
    python3-pip
COPY --from=prepare_requirements_txt_complete /app /app
RUN pip install -r /app/requirements.txt
COPY --from=prepare_build /app /app
WORKDIR /app
EXPOSE 3000
CMD ["npm", "run", "start"]
