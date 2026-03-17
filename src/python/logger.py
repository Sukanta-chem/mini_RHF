class Logger:

    def __init__(self, filename):
        self.file = open(filename, "w")

    def write(self, text):
        print(text)
        self.file.write(text + "\n")  

    def close(self):
        self.file.close()
