#-*- coding: utf-8 -*-
from __future__ import print_function
import sys
import xmltodict
import lxml.html
import unidecode
import json
from xml.parsers.expat import ExpatError
from importlib import reload

if sys.version_info[0] < 3:
    reload(sys)
    sys.setdefaultencoding('utf-8')

def getTargetitem(sourcedict, itemlist):
    #input hierarchically as (level1.level2.level3...)
    cur_item = sourcedict
    for i in itemlist:
        if(not isinstance(cur_item, dict)):
            break
        cur_item = cur_item.setdefault(i, "ERR")
    if(cur_item == None):
        return "ERR"
    return cur_item

def getPubdate(sourcedict):
    #Get publication dates from the data
    PubDate = getTargetitem(sourcedict, ['xocs:doc', 'xocs:item', 'item', 'bibrecord', 'head', 'source', 'publicationdate', 'year'])
    if(PubDate == "ERR"):
        PubDate = getTargetitem(sourcedict, ['xocs:doc', 'xocs:item', 'item', 'bibrecord', 'head', 'source', 'publicationdate', 'year'])
    if(PubDate == "ERR"): 
        PubDate = getTargetitem(sourcedict, ['xocs:doc', 'xocs:meta', 'xocs:sort-year'])     
    return PubDate

def getEID(sourcedict):
    EID = getTargetitem(sourcedict, ['xocs:doc', 'xocs:meta', 'xocs:eid'])
    return EID

def getSrcID(sourcedict):
    SrcID = getTargetitem(sourcedict, ['xocs:doc', 'xocs:item', 'item', 'bibrecord', 'head', 'source', '@srcid'])
    return SrcID

def getSrcTitle(sourcedict):
    SrcTitle = getTargetitem(sourcedict, ['xocs:doc', 'xocs:item', 'item', 'bibrecord', 'head', 'source', 'sourcetitle'])
    return SrcTitle

def getSrcType(sourcedict):
    SrcType = getTargetitem(sourcedict, ['xocs:doc', 'xocs:meta', 'xocs:srctype'])
    return SrcType

def getDocType(sourcedict):
    DocType = getTargetitem(sourcedict, ['xocs:doc', 'xocs:meta', 'cto:doctype'])
    return DocType

def getSrcISSN(sourcedict, isPrinted):
    ISSN = ""
    if(isPrinted):
        listISSN = getTargetitem(sourcedict, ['xocs:doc', 'xocs:item', 'item', 'bibrecord', 'head', 'source', 'issn'])
        if(isinstance(listISSN, list)):
            for vals in listISSN:
                if(vals.setdefault("@type", "###ERRMSG###") == "print"):
                    ISSN = vals.get("#text")
        elif(isinstance(listISSN, dict)):
            ISSN = listISSN.get("#text")
        else:
#            print(listISSN)
            pass 
         
    else:
        listISSN = getTargetitem(sourcedict, ['xocs:doc', 'xocs:item', 'item', 'bibrecord', 'head', 'source', 'issn'])
        if(isinstance(listISSN, list)):
            for vals in listISSN:
                if(vals.setdefault("@type", "###ERRMSG###") == "electronic"):
                    ISSN = vals.get("#text")
        elif(isinstance(listISSN, dict)):
            ISSN = listISSN.get("#text")
        else:
#            print(listISSN)
            pass
    if(ISSN == "ERR"):
        return ""
    else:
        return ISSN

def getASJC(sourcedict):
    allTempcode = getTargetitem(sourcedict, ['xocs:doc', 'xocs:item', 'item', 'bibrecord', 'head',
                                           'enhancement', 'classificationgroup', 'classifications'])
    if(allTempcode == "ERR"):
        ASJC = "-"
    elif(isinstance(allTempcode, list)):
        allCode = {}
        for i in allTempcode:
            allCode[i['@type']] = i['classification']
            ASJC = allCode.setdefault('ASJC', '-')
            if(isinstance(ASJC, list)):
                ASJC = ";".join(ASJC)
    else:
        Cur = allTempcode.setdefault('ASJC', '-')
        if(Cur == '-'):
            ASJC = "-"
        elif(isinstance(Cur, list)):
            ASJC = ";".join(Cur)
        else:
            ASJC = Cur
    return ASJC

def getMESH(sourcedict):
    allTempcode = getTargetitem(sourcedict, ['xocs:doc', 'xocs:item', 'item', 'bibrecord', 'head',
                                           'enhancement', 'descriptorgroup', 'descriptors'])
    if(allTempcode == "ERR"):
        MESH = "-"
    elif(isinstance(allTempcode, list)):
        allTerm = {}
        for i in allTempcode:
            allTerm[i['@type']] = i['descriptor']
        Cur = allTerm.setdefault('MSH', '-')
        if(Cur == '-'):
            MESH = "-"
        elif(isinstance(Cur, list)):
            MESH = ";".join([x['mainterm']['#text'] for x in Cur])
        else:
            MESH = Cur['mainterm']['#text']
    else:
        Cur = allTempcode.setdefault('MSH', '-')
        if(Cur == '-'):
            MESH = "-"
        elif(isinstance(Cur, list)):
            MESH = ";".join([x['mainterm']['#text'] for x in Cur])
        else:
            MESH = Cur['mainterm']['#text']

def getTitle(sourcedict):
    TitleList = getTargetitem(sourcedict, ['xocs:doc', 'xocs:item', 'item', 'bibrecord', 'head', 
                                          'citation-title', 'titletext', '#text'])
    if(isinstance(TitleList, unicode) or isinstance(TitleList, str)):
         return getStripTitle(TitleList)
    elif(isinstance(TitleList, list)):
        for val in TitleList:
            if(val.get('@original') == 'y'):
                return getStripTitle(val['#text'])
        for val in TitleList:
            if(val.get('xml:lang') == 'eng'):
                return getStripTitle(val['#text'])
    else:
        return ""

def getRefEID(sourcedict):
    RefList = getTargetitem(sourcedict, ['xocs:doc', 'xocs:meta', 'cto:ref-id'])
    if(RefList == "ERR"):
        return -1
    if(not isinstance(RefList, list)):
        RefList = [RefList]
    return RefList

def getStripTitle(orgtitle):
    nowtitle = unidecode.unidecode(unicode(orgtitle.lower().strip()))
    filteredtitle = "".join([x for x in nowtitle if unicode(x).isalpha() or unicode(x).isnumeric() or unicode(x).isspace()])
    return filteredtitle

def getRegauthorkey(AuthorKeys):
    joinAuthorKey = ""
    NewList = []
    if(isinstance(AuthorKeys, list)):
        for i in range(len(AuthorKeys)):
            if(isinstance(AuthorKeys[i], dict)):
                if("#text" in AuthorKeys[i]):
                    NewList.append(getRegtext(AuthorKeys[i]['#text']))
                else:
                    pass
            else:
                NewList.append(AuthorKeys[i])
        joinAuthorKey = ":::".join([x for x in NewList if x is not None])
    elif(isinstance(AuthorKeys, str)):
        if(AuthorKeys=="ERR"):
            joinAuthorKey == ""
        else:
            joinAuthorKey == getRegtext(AuthorKeys)
    if(joinAuthorKey is not None):
        return joinAuthorKey.encode('utf8')
    else:
        return ""

def getRegtext(text):
    return " ".join("".join(text.splitlines()).split())

def getDOI(sourcedict):
    DOI = getTargetitem(sourcedict, ['xocs:doc', 'xocs:meta', 'xocs:doi'])
    if(DOI=="ERR"):
        DOI = getTargetitem(sourcedict, ['xocs:doc', 'xocs:item', 'item', 'bibrecord', 'item-info', 'itemidlist', 'ce:doi'])
    if(DOI=="ERR"):
        DOI = ""
    return DOI

def getAuthorInfoJSON(sourcedict):
    author_group = getTargetitem(sourcedict, ['xocs:doc', 'xocs:item', 'item', 'bibrecord', 'head', 'author-group'])
    infors = []
    #print(author_group)
    if(not isinstance(author_group, list)):
        infors.append(json.dumps(author_group))
    else:
        for group in author_group :
            infors.append(json.dumps(group))
    return ";".join(infors)

def getAffsAuthorinfo(sourcedict):
    # Return the information of Affilations with author ID
    author_group = getTargetitem(sourcedict, ['xocs:doc', 'xocs:item', 'item', 'bibrecord', 'head', 'author-group'])
    af_info = []
    if(not isinstance(author_group, list)):
        cur_aff_info = {}
        author_info = []
        ctry = ""
        afid = ""
        if(author_group == "ERR"):
            pass
        elif('affiliation' not in author_group):
            pass
        else:
            if('@country' not in author_group['affiliation']):
                pass
            else:
                ctry = author_group['affiliation']['@country']
            if('@afid' not in author_group['affiliation']):
                pass
            else:
                afid = author_group['affiliation']['@afid']
        
        #Get Author IDS for the affiliations
        if(author_group == "ERR"):
            pass
        else:
            author_list = author_group['author']
            if(isinstance(author_list, list)):
                for aut in author_list:
                    cur_author_auid = ""
                    cur_author_indexed_name = ""
                    cur_author_surname = ""
                    cur_author_given_name = ""
                    cur_author_sequence = ""
                    cur_author_info = {}
                    if('@auid' not in aut):
                        pass
                    else:
                        cur_author_auid = (aut['@auid'])
                    if('@seq' not in aut):
                        pass
                    else:
                        cur_author_sequence = (aut['@seq'])                         
                    if('ce:indexed-name' not in aut):
                        pass
                    else:
                        cur_author_indexed_name = (aut['ce:indexed-name'])                             
                    if('ce:surname' not in aut):
                        pass
                    else:
                        cur_author_surname = (aut['ce:surname'])    
                    if('ce:given-name' not in aut):
                        pass
                    else:
                        cur_author_given_name = (aut['ce:given-name'])    
                    cur_author_info["auid"] = cur_author_auid
                    cur_author_info["indexedname"] = cur_author_indexed_name
                    cur_author_info["surname"] = cur_author_surname
                    cur_author_info["givenname"] = cur_author_given_name
                    cur_author_info["seq"] = cur_author_sequence
                    author_info.append(cur_author_info)

            else:
                cur_author_auid = ""
                cur_author_indexed_name = ""
                cur_author_surname = ""
                cur_author_given_name = ""
                cur_author_sequence = ""
                cur_author_info = {}                
                if('@auid' not in author_list):
                    pass
                else:
                    cur_author_auid = (author_list['@auid'])
                if('@seq' not in author_list):
                    pass
                else:
                    cur_author_sequence = (author_list['@seq'])                    
                if('ce:indexed-name' not in author_list):
                    pass
                else:
                    cur_author_indexed_name = (author_list['ce:indexed-name'])    
                if('ce:surname' not in author_list):
                    pass
                else:
                    cur_author_surname = (author_list['ce:surname'])    
                if('ce:given-name' not in author_list):
                    pass
                else:
                    cur_author_given_name = (author_list['ce:given-name'])
                cur_author_info["auid"] = cur_author_auid
                cur_author_info["indexedname"] = cur_author_indexed_name
                cur_author_info["surname"] = cur_author_surname
                cur_author_info["givenname"] = cur_author_given_name
                cur_author_info["seq"] = cur_author_sequence
                author_info.append(cur_author_info)                       

        cur_aff_info["afid"] = afid
        cur_aff_info["ctry"] = ctry       
        cur_aff_info["authors"] = author_info
        af_info.append(cur_aff_info)
        
    else:
        for group in author_group :
            cur_aff_info = {}
            author_info = []
            ctry = ""
            afid = ""
            if(not isinstance(group, dict)):
                continue
            elif('affiliation' not in group):
                pass
            else:
                if('@country' not in group['affiliation']):
                    pass
                else:
                    ctry = group['affiliation']['@country']
                if('@afid' not in group['affiliation']):
                    pass
                else:
                    afid = group['affiliation']['@afid']
            
            #Get Author IDS for the affiliations
            if(not isinstance(group, dict)):
                continue
            else:
                if('author' not in group):
                    continue
                author_list = group['author']
                if(isinstance(author_list, list)):
                    cur_author_auid = ""
                    cur_author_indexed_name = ""
                    cur_author_surname = ""
                    cur_author_given_name = ""
                    cur_author_sequence = ""
                    cur_author_info = {}
                    for aut in author_list:
                        if('@auid' not in aut):
                            pass
                        else:
                            cur_author_auid = (aut['@auid'])
                        if('@seq' not in aut):
                            pass
                        else:
                            cur_author_sequence = (aut['@seq'])                              
                        if('ce:indexed-name' not in aut):
                            pass
                        else:
                            cur_author_indexed_name = (aut['ce:indexed-name'])    
                        if('ce:surname' not in aut):
                            pass
                        else:
                            cur_author_surname = (aut['ce:surname'])    
                        if('ce:given-name' not in aut):
                            pass
                        else:
                            cur_author_given_name = (aut['ce:given-name'])    
                        cur_author_info["auid"] = cur_author_auid
                        cur_author_info["indexedname"] = cur_author_indexed_name
                        cur_author_info["surname"] = cur_author_surname
                        cur_author_info["givenname"] = cur_author_given_name
                        cur_author_info["seq"] = cur_author_sequence
                        author_info.append(cur_author_info)

                else:
                    cur_author_auid = ""
                    cur_author_indexed_name = ""
                    cur_author_surname = ""
                    cur_author_given_name = ""
                    cur_author_sequence = ""                    
                    cur_author_info = {}                
                    if('@auid' not in author_list):
                        pass
                    else:
                        cur_author_auid = (author_list['@auid'])
                    if('@seq' not in author_list):
                        pass
                    else:
                        cur_author_sequence = (author_list['@seq'])                         
                    if('ce:indexed-name' not in author_list):
                        pass
                    else:
                        cur_author_indexed_name = (author_list['ce:indexed-name'])    
                    if('ce:surname' not in author_list):
                        pass
                    else:
                        cur_author_surname = (author_list['ce:surname'])    
                    if('ce:given-name' not in author_list):
                        pass
                    else:
                        cur_author_given_name = (author_list['ce:given-name'])
                    cur_author_info["auid"] = cur_author_auid
                    cur_author_info["indexedname"] = cur_author_indexed_name
                    cur_author_info["surname"] = cur_author_surname
                    cur_author_info["givenname"] = cur_author_given_name
                    cur_author_info["seq"] = cur_author_sequence
                    author_info.append(cur_author_info)        
                    
            cur_aff_info["afid"] = afid
            cur_aff_info["ctry"] = ctry           
            cur_aff_info["authors"] = author_info
            af_info.append(cur_aff_info)
    #print(json.dumps(af_info))
    
    return json.dumps(af_info)
            
def getAuid(sourcedict):
    author_group = getTargetitem(sourcedict, ['xocs:doc', 'xocs:item', 'item', 'bibrecord', 'head', 'author-group'])
    auids = []
    if(not isinstance(author_group, list)):
        if(author_group == "ERR"):
            return ""
        if('author' not in author_group):
            return ""
        author_list = author_group['author']
        if(isinstance(author_list, list)):
            for aut in author_list:
                if('@auid' not in aut):
                    continue
                auids.append(aut['@auid'])
        else:
            if('@auid' not in author_group['author']):
                return ""
            auids.append(author_group['author']['@auid'])
    else:
        for group in author_group :
            if(not isinstance(group, dict)):
                continue
            if('author' not in group):
                continue
            if(isinstance(group['author'], list)):
                for author in group['author'] :
                    if('@auid' not in author):
                        continue
                    auids.append(author['@auid'])
            else:
                if('@auid' not in group['author']):
                    continue
                auids.append(group['author']['@auid'])
    return ";".join(auids)

def getAfid(sourcedict):
    author_group = getTargetitem(sourcedict, ['xocs:doc', 'xocs:item', 'item', 'bibrecord', 'head', 'author-group'])
    afids = []
    if(not isinstance(author_group, list)):
        if(author_group == "ERR"):
            return ""
        if('affiliation' not in author_group):
            return ""
        if('author' not in author_group):
            return ""            
        for i in range(len(author_group['author'])):
            if('@afid' not in author_group['affiliation']):
                return ""
            afids.append(author_group['affiliation']['@afid'])
    else:
        for group in author_group :
            if(not isinstance(group, dict)):
                continue        
            if('affiliation' not in group):
                continue
            if('@afid' not in group['affiliation']):
                continue
            if('author' not in group):
                continue
            for i in range(len(group['author'])):
                if('@auid' not in group['author']):
                    continue            
                afids.append(group['affiliation']['@afid'])
    return ";".join(afids)

def getAuthorcountry(sourcedict):
    author_group = getTargetitem(sourcedict, ['xocs:doc', 'xocs:item', 'item', 'bibrecord', 'head', 'author-group'])
    auct = []
    if(not isinstance(author_group, list)):
        if(author_group == "ERR"):
            return ""
        if('affiliation' not in author_group):
            return ""
        if('author' not in author_group):
            return ""              
        for i in range(len(author_group['author'])):
            if('@country' not in author_group['affiliation']):
                return ""
            auct.append(author_group['affiliation']['@country'])
    else:
        for group in author_group :
            if(not isinstance(group, dict)):
                continue        
            if('affiliation' not in group):
                continue
            if('@country' not in group['affiliation']):
                continue
            if('author' not in group):
                continue
            for i in range(len(group['author'])):
                auct.append(group['affiliation']['@country'])
    return ";".join(auct)
    
def openXml(xmlstring):
    try:
        tempdict = xmltodict.parse(xmlstring)
    except ExpatError as e:
        xmlstring = fixBadxml(xmlstring)
        tempdict = xmltodict.parse(xmlstring)
    return tempdict

def fixBadxml(xmlstring):
    try:
        lx = lxml.html.fromstring(xmlstring)
    except Exception as e:
        return "<xocs:doc><error>badxml</error></xocs:doc>"
    fixed = lxml.html.tostring(lx, method='xml')
    return fixed

if(__name__ == "__main__"):
    xmlbuffer = ""
    count = 0
    print("\t".join(["#eid", "pubyear", "asjc", "srcid", "srctitle", "srctype", "doctype", "authorinfo"]))
    for line in sys.stdin:
        lines = line.split("<?xml")
        xmlbuffer += lines[0]
        if(line.find("</xocs:doc") >= 0):
            tempdict = openXml(xmlbuffer)
            EID = getEID(tempdict)
            PubYear = getPubdate(tempdict)
            #ASJC = getASJC(tempdict)
            #SrcID = getSrcID(tempdict)
            #SrcTitle = getSrcTitle(tempdict)
            #SrcType = getSrcType(tempdict)
            #DocType = getDocType(tempdict)
            #AuthorINFO = getAuthorInfoJSON(tempdict)
            #AuthorINFO = getAffsAuthorinfo(tempdict)
            AuthorCountry = getAuthorcountry(tempdict)
            #print("Parsed")
            print(EID, PubYear, AuthorCountry, sep="\t")
            xmlbuffer = ""