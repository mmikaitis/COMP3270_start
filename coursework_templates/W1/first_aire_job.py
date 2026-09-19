#!/usr/bin/env python3

import os
import socket
from datetime import datetime

def main():
    # Get username
    user = os.environ.get("USER")

    # Get node hostname
    node = socket.gethostname()

    # Get IP of the node
    ip = socket.gethostbyname(socket.gethostname())

    # Get current date and time
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    print(f"User who submitted this job: {user}")
    print(f"Node that this job run on: {node}")
    print(f"IP: {ip}")
    print(f"Date & Time: {now}")

if __name__ == "__main__":
    main()

