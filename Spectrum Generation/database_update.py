#!/usr/bin/python
# -*- coding: utf-8 -*-
import xdrlib, sys
import os
import xlrd
import queue
import urllib
import urllib.request as urllib2
import re
import threading
import time


def download_from_hmbd(name):
    # name = 'Nicotinic acid'
    name_quote = urllib.parse.quote_plus(name)
    # print(name)
    url = 'http://www.hmdb.ca/unearth/q?utf8=%E2%9C%93&query=' + name_quote + '&searcher=metabolites&button='
    # print(url)
    content = urllib2.urlopen(url).read().decode()

    m = re.search("hit-name.{40}", content)
    m = re.split('"', m.group())
    ID = re.split('/', m[2])[-1]

    url = 'http://www.hmdb.ca/metabolites/' + ID
    content = urllib2.urlopen(url).read().decode()
    m = re.search(".{40}1H NMR Spectrum", content)
    try:
        nmrID = re.search("/spectra/nmr_one_d/[0-9]*", m.group()).group()
    except:
        print("\nCan not find the followed metabolite")
        print(name)
            


    url = 'http://www.hmdb.ca' + nmrID
    content = urllib2.urlopen(url).read().decode()
    m = re.search("List of chemical shift values for the spectrum.{200}", content)
    m = re.search("a href=.*", m.group())
    file_url = re.split('"', m.group())[1]

    data = urllib2.urlopen(file_url).read().decode()

    f = open('database/' + name + '.txt', 'w')
    f.write(data)
    f.close()
    # print(content)


def open_excel(file='metabolites.xlsx'):
    try:
        data = xlrd.open_workbook(file)
        return data
    except Exception as e:
        print(str(e))


# 根据索引获取Excel表格中的数据   参数:file：Excel文件路径     colnameindex：表头列名所在行的所以  ，by_index：表的索引
def excel_table_byindex(file='metabolites.xlsx', colnameindex=0, by_index=0):
    data = open_excel(file)
    table = data.sheets()[by_index]
    nrows = table.nrows  # 行数
    ncols = table.ncols  # 列数
    colnames = table.row_values(colnameindex)  # 某一行数据
    list = []
    for rownum in range(1, nrows):

        row = table.row_values(rownum)
        if row:
            app = {}
            for i in range(len(colnames)):
                app[colnames[i]] = row[i]
            list.append(app)
    return list

class myThread (threading.Thread):   #继承父类threading.Thread
    def __init__(self, threadID, name, counter):
        threading.Thread.__init__(self)
        self.threadID = threadID
        self.name = name
        self.counter = counter
    def run(self):                   #把要执行的代码写到run函数里面 线程在创建后会直接运行run函数 
        print( "\nStarting " + self.name)
        download_from_hmbd(self.name)
        print( "\nExiting " + self.name)

def main():
    tables = excel_table_byindex()
    q = queue.Queue()
    for row in tables:
        name_check = row['name']
        name_path = os.getcwd()+ '\\database\\'+ name_check+'.txt'
        print('Check ' +name_check + '...')
        if(os.path.exists(name_path)):
            continue
        else:
            q.put(name_check )

        #print(row)
        #print(row['name'] )
        
    while not q.empty():
        name_uodate = q.get()
        thread1 = myThread(1, name_uodate, 1).start()
        #download_from_hmbd(name_check)

    print( "Exiting Main Thread")


if __name__ == "__main__":
    main()
    print("Update finish!")
