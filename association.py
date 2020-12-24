
import numpy as np
import math
import time
import sys


#input txt file name
#output a Matrix called sequence_l where each row represent sequence of gene (sequence)
def preprocess_2(file):
  f = open(file).readlines()
  sequence_list = []
  for lines in f:
    line = lines.split("\t")
    gene_l = []
    for x in range (len(line)-1):
      gene_l.append("G"+str(x+1)+"_"+line[x])
    gene_l.append(line[-1].strip("\n"))
    sequence_list.append(gene_l)
  return sequence_list

#given a sequence and a list of sequence
#output number of appearence of the sequence within the list
def get_support_count2(sequence,sequence_l):
  sequence_set = set(sequence.split(','))
  count = 0
  for s in sequence_l:
    b = True
    for y in sequence_set:
      y_index = y[1:y.index('_')]
      index = int(y_index)
      if (s[index-1]!=y):
        b = False
        break
    if (b):
      count = count +1
  return count

#generate candidate set
#input gene martix, size of candidate set desire, previous k-1 size frequenset, if size = 0 then feed with empty list 
def generate_set(sequence_l,size,previous_frequent_set):
  candidate_set = []
  if (size == 1):
    for sequence in sequence_l:
      for x in range (len(sequence)-1):
        if (sequence[x] not in candidate_set):
          candidate_set.append(sequence[x])
    return [elem for elem in candidate_set]
  else:
    for x in range (len(previous_frequent_set)-1) :
      for y in range (x+1,len(previous_frequent_set)):
        #print('prev ',previous_frequent_set[x])
        set1 = set(previous_frequent_set[x].split(','))
        set2 = set(previous_frequent_set[y].split(','))
        if (len(set1.union(set2))==size):
          l = list(set1.union(set2))
          l.sort()
          s = ",".join(l)
          #if (s not in candidate_set):
            #print('s',s)
          candidate_set.append(s)
  return list(set(candidate_set))

#input gene matrix, candidate sets, minimum support threshhold, previous infrequent_set
def scan_and_generate_frequenset(sequence_l,candidates,min_sup_count,prev_infrequent_set):
  freq_set = []
  infrequent_set  = []
  count_set = []
  #pruning
  count = 0
  for candidate in candidates:
    # for item in prev_infrequent_set:
    #   l = item.split(",")
    #   if(all(x in candidate for x in l)): 
    #     infrequent_set.append(candidate)
    #     continue
    count = get_support_count2(candidate,sequence_l)
    if (count>=min_sup_count):
      freq_set.append(candidate)
      count_set.append(count)
    else:
      infrequent_set.append(candidate)
  return freq_set,prev_infrequent_set+infrequent_set,count_set


#this function generates the frequent sets of all length
#input min support thresh hold and gene matrix
#out frequenset and its support count
def Apriori (min_sup,sequence_l):
  start = time.time()
  b = False
  min_sup_count = math.ceil(min_sup * len(sequence_l))
  size =1
  frequent_set =[]
  inf_set = []
  total = 0
  ret_freq_set = []
  ret_count_set =[]
  while (not b):
    candidate_set = generate_set(sequence_l,size,frequent_set)
    frequent_set,inf_set,count_set = scan_and_generate_frequenset(sequence_l,candidate_set,min_sup_count,[])
    ret_freq_set = ret_freq_set + frequent_set
    ret_count_set =ret_count_set+count_set
    if(len(frequent_set) == 0):
      b = False
      break
    print("number of length-",size," frequent itemsets: " ,len(frequent_set))
    size = size +1
    total = total + len(frequent_set)
  print("number of all lengths frequent itemsets: " ,total)
  return ret_freq_set , ret_count_set



def rule_on_one_freq_set(head,body,sup_count,min_confidence,sequence_l):
    head_l = head.split(',')
    size = len(head_l)
    #case size one frequent set
    if (size==1 and body == ""):
      return [],[]
    head_count = get_support_count2(head,sequence_l)
    confidence = sup_count/head_count
    ret_head_l = []
    ret_body_l = []
    if (confidence<min_confidence):
      return [],[]
    elif (size == 1):
      ret_head_l.append(head)
      ret_body_l.append(body)
      return ret_head_l,ret_body_l
    else:
      if (body !=""):
        ret_head_l.append(head)
        ret_body_l.append(body)
      #recursively generate all sub rule
      for i in range (size):
        #form new body and head
        sub_head_l = head_l[0:i] +head_l[i+1:]
        if (body == ""):
          body_l = head_l[i]
        else:
          l = body +","+head_l[i]
          l = l.split(',')
          l.sort()
          body_l =  ",".join(l)
        #recursive call
        sub_rule_head, sub_rule_body = rule_on_one_freq_set(",".join(sub_head_l), body_l ,sup_count,min_confidence,sequence_l)
        ret_head_l = ret_head_l +sub_rule_head
        ret_body_l = ret_body_l + sub_rule_body
      ret_head_l_final = []
      ret_body_l_final = []
      #delete duplicates
      for i in range (len(ret_head_l)):
        if (ret_head_l[i] not in ret_head_l_final):
          ret_head_l_final.append(ret_head_l[i])
          ret_body_l_final.append(ret_body_l[i])
      return ret_head_l_final,ret_body_l_final

def rule_generation(min_sup, min_confidence,sequence_l):
  freq_set, count_set = Apriori(min_sup,sequence_l)
  head = []
  body = []
  for i in range (len(freq_set)):
    h,b = rule_on_one_freq_set(freq_set[i],'',count_set[i],min_confidence,sequence_l)
    head = head + h
    body = body + b
  return head ,body

def template1(arg1,arg2,arg3,min_sup, min_confidence,sequence_l):
  count = 0
  head,body = rule_generation(min_sup, min_confidence,sequence_l)
  h = []
  b = []
  if (arg1=="RULE"):
    if (arg2=="ANY"):
      for x in range (len(head)):
        if (any(item in head[x] for item in arg3)):
          count = count + 1
          h.append(head[x])
          b.append(body[x])
        elif (any(item in body[x] for item in arg3)):
          count = count + 1
          h.append(head[x])
          b.append(body[x])
    elif (arg2 =="NONE"):
      for x in range (len(head)):
        if (all(item not in head[x] for item in arg3) and all(item not in body[x] for item in arg3)):
          count = count + 1
          h.append(head[x])
          b.append(body[x])
    else:
      arg2 = int(arg2)
      for x in range (len(head)):
        req_count = 0
        for item in arg3:
          if (item in head[x]):
            req_count = req_count +1
          elif (item in body[x]):
            req_count = req_count +1
        if (req_count == arg2):
          count = count + 1
          h.append(head[x])
          b.append(body[x])
  elif (arg1=="HEAD"):
    if (arg2=="ANY"):
      for x in range (len(head)):
        if (any(item in head[x] for item in arg3)):
          count = count + 1
          h.append(head[x])
          b.append(body[x])
    elif (arg2 =="NONE"):
      for x in range (len(head)):
        if (all(item not in head[x] for item in arg3)):
          count = count + 1
          h.append(head[x])
          b.append(body[x])
    else:
      arg2 = int(arg2)
      for x in range (len(head)):
        req_count = 0
        for item in arg3:
          if (item in head[x]):
            req_count = req_count +1
        if (req_count == arg2):
          count = count + 1
          h.append(head[x])
          b.append(body[x])
  elif (arg1=="BODY"):
    if (arg2=="ANY"):
      for x in range (len(body)):
        if (any(item in body[x] for item in arg3)):
          count = count + 1
          h.append(head[x])
          b.append(body[x])
    elif (arg2 =="NONE"):
      for x in range (len(body)):
        if (all(item not in body[x] for item in arg3)):
          count = count + 1
          h.append(head[x])
          b.append(body[x])
    else:
      arg2 = int(arg2)
      for x in range (len(body)):
        req_count = 0
        for item in arg3:
          if (item in body[x]):
            req_count = req_count +1
        if (req_count == arg2):
          count = count + 1
          h.append(head[x])
          b.append(body[x])
  else:
    print("invalid input")
  return h,b,count

def template2(arg1,arg2,min_sup, min_confidence,sequence_l):
  count = 0
  h = []
  b = []
  arg2 = int(arg2)
  head,body = rule_generation(min_sup, min_confidence,sequence_l)
  if (arg1 == "RULE"):
    for x in range(len(head)):
      s = head[x] + body[x]
      size = s.count('_')
      if (size >= arg2):
        count = count +1
        h.append(head[x])
        b.append(body[x])
  elif (arg1 == "HEAD"):
    for x in range(len(head)):
      s = head[x] 
      size = s.count('_')
      if (size >= arg2):
        count = count +1
        h.append(head[x])
        b.append(body[x])
  elif (arg1 == "BODY"):
    for x in range(len(body)):
      s = body[x]
      size = s.count('_')
      if (size >= arg2):
        count = count +1
        h.append(head[x])
        b.append(body[x])
  return h,b,count

def template3(logical_op,arg1_1,arg1_2,arg1_3,arg2_1,arg2_2,arg2_3,min_sup, min_confidence,sequence_l):
  if (logical_op== "1or1"):
    h1,b1,c1 = template1(arg1_1,arg1_2,arg1_3,min_sup, min_confidence,sequence_l)
    h = h1
    b = b1
    c = c1
    r = []
    for i in range (len(h1)):
      r.append(h1[i]+b1[i])
    h2,b2,c2 = template1(arg2_1,arg2_2,arg2_3,min_sup, min_confidence,sequence_l)
    for x in range (len(h2)):
      if (h2[x]+b2[x] not in r):
        h.append(h2[x])
        b.append(b2[x])
        c = c+1
  elif (logical_op== "1and1"):
    h = []
    b = []
    c = 0
    h1,b1,c1 = template1(arg1_1,arg1_2,arg1_3,min_sup, min_confidence,sequence_l)
    r = []
    for i in range (len(h1)):
      r.append(h1[i]+b1[i])
    h2,b2,c2 = template1(arg2_1,arg2_2,arg2_3,min_sup, min_confidence,sequence_l)
    for x in range (len(h2)):
      if (h2[x]+b2[x] in r):
        h.append(h2[x])
        b.append(b2[x])
        c = c+1
  elif (logical_op== "1or2"):
    print (arg1_1,arg1_2,arg1_3,min_sup, min_confidence)
    print(arg2_1,arg2_2)
    h1,b1,c1 = template1(arg1_1,arg1_2,arg1_3,min_sup, min_confidence,sequence_l)
    h = h1
    b = b1
    c = c1
    r = []
    for i in range (len(h1)):
      r.append(h1[i]+b1[i])
    h2,b2,c2 = template2(arg2_1,arg2_2,min_sup, min_confidence,sequence_l)
    for x in range (len(h2)):
      if (h2[x]+b2[x] not in r):
        h.append(h2[x])
        b.append(b2[x])
        c = c+1
  elif (logical_op== "1and2"):
    h = []
    b = []
    c = 0
    h1,b1,c1 = template1(arg1_1,arg1_2,arg1_3,min_sup, min_confidence,sequence_l)
    h2,b2,c2 = template2(arg2_1,arg2_2,min_sup, min_confidence,sequence_l)
    r = []
    for i in range (len(h1)):
      r.append(h1[i]+b1[i])
    for x in range (len(h2)):
      if (h2[x]+b2[x] in r):
        h.append(h2[x])
        b.append(b2[x])
        c = c+1
  elif (logical_op== "2or2"):
    h1,b1,c1 = template2(arg1_1,arg1_2,min_sup, min_confidence,sequence_l)
    h = h1
    b = b1
    c = c1
    r = []
    for i in range (len(h1)):
      r.append(h1[i]+b1[i])
    h2,b2,c2 = template2(arg2_1,arg2_2,min_sup, min_confidence,sequence_l)
    for x in range (len(h2)):
      if (h2[x]+b2[x] not in r):
        h.append(h2[x])
        b.append(b2[x])
        c = c+1
  elif (logical_op== "2and2"):
    h = []
    b = []
    c = 0
    h1,b1,c1 = template2(arg1_1,arg1_2,min_sup, min_confidence,sequence_l)
    r = []
    for i in range (len(h1)):
      r.append(h1[i]+b1[i])
    h2,b2,c2 = template2(arg2_1,arg2_2,min_sup, min_confidence,sequence_l)
    for x in range (len(h2)):
      if (h2[x]+b2[x] in r):
        h.append(h2[x])
        b.append(b2[x])
        c = c+1
  else: 
    print('invalid operator')
    return None
  return h,b,c

start = time.time()
#f = 'association-rule-test-data.txt'
f = sys.argv[1]
min_sup = float(sys.argv[2])
sequence_l = preprocess_2(f)

if (len(sys.argv)==3):
  Apriori(min_sup,sequence_l)
else:
  min_confidence = float(sys.argv[3])
  if (sys.argv[4]=="template1"):
    arg1 = sys.argv[5]
    arg2 = sys.argv[6]
    arg3 = sys.argv[7].split(',')
    h,b,count = template1(arg1,arg2,arg3,min_sup, min_confidence,sequence_l)
    print('result,count = asso_rule.template1(',arg1,',',arg2,',',arg3,'): ',count)
  elif(sys.argv[4]=="template2"):
    arg1 = sys.argv[5]
    arg2 = sys.argv[6]
    h,b,count = template2(arg1,arg2,min_sup, min_confidence,sequence_l)
    print('result,count = asso_rule.template2(',arg1,',',arg2,'): ',count)
  elif(sys.argv[4]=="template3"):
    logical_op = sys.argv[5]
    arg1_1 = sys.argv[6]
    arg1_2 = sys.argv[7]
    arg1_3 = sys.argv[8].split(',')
    arg2_1 = sys.argv[9]
    arg2_2 = sys.argv[10]
    arg2_3 = sys.argv[11].split(',')
    h,b,count = template3(logical_op,arg1_1,arg1_2,arg1_3,arg2_1,arg2_2,arg2_3,min_sup, min_confidence,sequence_l)
    print('result,count = asso_rule.template3(',logical_op,',', arg1_1,',',arg1_2,',',arg1_3, ',',arg2_1,',',arg2_2,',',arg2_3,'): ',count)


def print_rule(h,b):
  for x in range (len(h)):
    print('head',h[x], ' body: ',b[x])




