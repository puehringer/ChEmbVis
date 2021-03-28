FROM node:12-alpine
RUN mkdir /frontend
WORKDIR /frontend

# Install app dependencies
COPY package.json ./
# COPY package-lock.json ./
RUN npm install

# Add app
COPY . ./

# Run the dev server (or install deps if node_modules is not existent)
CMD [ -d "node_modules" ] && npm run start || npm i && npm start
