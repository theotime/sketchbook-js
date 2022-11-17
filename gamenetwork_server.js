// npm init
// npm install express@4
// npm install socket.io

// node index.js

/*
    Simple relay server for a peer to peer network.
*/

const express = require("express");
const app = express();

const http = require("http");
const server = http.createServer(app);

const { Server } = require("socket.io");
const io = new Server(server);

const port = 5000;

app.get("/", (req, res) => {
    res.sendFile(__dirname + "/index.html");
});


// AUTH
io.on("connection", (socket) => {

    console.log("user connected.");
    socket.emit("auth-ack", socket.id);

    // RECV
    let serializeState = (id, state) => id + ":" + state;
    socket.on("update", (msg) => {
        //console.log("new state received.");
        io.emit("broadcast-update", serializeState(socket.id, msg));
    });

    socket.on("disconnect", () => {
        console.log("user disconnected.");
        io.emit("broadcast-disconnection", socket.id);
    });

});

// LISTEN
server.listen(port, () => {
    console.log("listening.");
});