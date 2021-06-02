from Bio import SeqIO
import sys,os,re
from collections import OrderedDict

def check_inter(data_list):
    flag = 0
    for i in range(len(data_list)):
        if data_list != []:
            key1 = data_list.pop()
            start1 = int(key1.split('|')[1])
            end1 = int(key1.split('|')[2])
            strand1 = key1.split('|')[3]
            for key in data_list:
                start2 = int(key.split('|')[1])
                end2 = int(key.split('|')[2])
                strand2 = key.split('|')[3]
                if start1 < end2 and end1 > start2 and strand1 == strand2:
                    flag = 1
                    break
    return flag
    
def multi_merge(new_data_dict):
    second_key_list = list(new_data_dict.keys())
    new_data_dict2 = new_data_dict
    new_data_dict2 = {}
    temp_list = []
    for i in range(len(second_key_list)):
        if second_key_list != []:
            key1 = second_key_list.pop()
            name = key1.split('|')[0]
            start1 = int(key1.split('|')[1])
            end1 = int(key1.split('|')[2])
            strand = key1.split('|')[3]
            for key in second_key_list:
                start2 = int(key.split('|')[1])
                end2 = int(key.split('|')[2])
                strand2 = key.split('|')[3]
                if start1 < end2 and end1 > start2 and strand == strand2:
                    if start1 <= start2:
                        start = start1
                    else:
                        start = start2
                    if end1 >= end2:
                        end = end1
                    else:
                        end = end2
                    new_key = name + '|' + str(start) + '|' + str(end) + '|' + strand
                    new_info = new_data_dict[key] + new_data_dict[key1]
                    new_data_dict2[new_key] = new_info
                    if key not in temp_list:
                        temp_list.append(key)
                    temp_list.append(key1)
                    second_key_list.remove(key)
        else:
            break
    for k,v in new_data_dict.items():
        if k not in temp_list:
            new_data_dict2[k] = v
    #调整字典的键
    new_data_dict3 = OrderedDict()
    for data,info in new_data_dict2.items():
        min = 999999999
        max = 0
        for each_data in info:
            start = int(each_data.split('|')[1])
            end = int(each_data.split('|')[2])
            strand0 = each_data.split('|')[3]
            if start <= min:
                min = start
            if end >= max:
                max = end
        new_data_key = data.split('|')[0] + '|' + str(min) + '|' + str(max) + '|' + strand0
        new_data_dict3[new_data_key] = info
    return new_data_dict3

def find_overlap(outdir,feature_type,source):
    data_list = []
    input = outdir+'/tmp.gb'
    if os.path.exists(input):
        #从genbank文件构建字典
        for record1 in SeqIO.parse(input, "genbank"):
            for feature in record1.features:
                strand = str(feature.strand)
                if feature.type == feature_type:
                    if feature.qualifiers['source'][0] == source:
                        name = str(feature.qualifiers['gene'][0])
                        pattern = "\d+"
                        result = re.findall(pattern,str(feature.location))
                        start = str(int(result[0]) + 1)
                        end = result[1]
                        start,end = exchange(start,end,strand)
                        pos = start + ":" + end
                        data = name+'|'+start+'|'+end+'|'+strand
                        data_list.append(data)
            break
        if data_list:
            #根据overlap分组
            data_dict = OrderedDict()
            while 1:
                if data_list != []:
                    data = data_list.pop()
                    start = int(data.split('|')[1])
                    end = int(data.split('|')[2])
                    strand1 = data.split('|')[3]
                    data_dict[data] = []
                    data_dict[data].append(data)
                    remove_list = []
                    for j in range(len(data_list)):
                        start_j = int(data_list[j].split('|')[1])
                        end_j = int(data_list[j].split('|')[2])
                        strand_j = data_list[j].split('|')[3]
                        if start < end_j and end > start_j and strand_j == strand1:
                            data_dict[data].append(data_list[j])
                            remove_list.append(data_list[j])
                    for li in remove_list:
                        data_list.remove(li)
                else:
                    break
            #调整字典的键
            new_data_dict = OrderedDict()
            for data,info in data_dict.items():
                min = 99999999999
                max = 0
                for each_data in info:
                    start = int(each_data.split('|')[1])
                    end = int(each_data.split('|')[2])
                    strand0 = each_data.split('|')[3]
                    if start <= min:
                        min = start
                    if end >= max:
                        max = end
                new_data_key = data.split('|')[0] + '|' + str(min) + '|' + str(max) + '|' + strand0
                new_data_dict[new_data_key] = info
            #多次合并
            new_data_dict1 = multi_merge(new_data_dict)
            flag = check_inter(list(new_data_dict1.keys()))
            if flag == 1:
                new_data_dict2_1 = multi_merge(new_data_dict1)
                flag1 = check_inter(list(new_data_dict2_1.keys()))
                if flag1 == 1:
                    new_data_dict2_2 = multi_merge(new_data_dict2_1)
                    flag2 = check_inter(list(new_data_dict2_2.keys()))
                    if flag2 == 1:
                        new_data_dict2 = multi_merge(new_data_dict2_2)
                    else:
                        new_data_dict2 = new_data_dict2_2
                else:
                    new_data_dict2 = new_data_dict2_1
            else:
                new_data_dict2 = new_data_dict1
        else:
            new_data_dict2 = {}
        return new_data_dict2
                
def get_data(outdir):
    resistance_dict,IS_dict,integron_dict,transposon_dict = OrderedDict(),OrderedDict(),OrderedDict(),OrderedDict()
    input = outdir+'/tmp.gb'

    if os.path.exists(input):

        for rec in SeqIO.parse(input,'genbank'):
            seq = str(rec.seq)
            endpos = str(len(seq))
        IS_dict = find_overlap(outdir,'mobile_element',"IS")
        integron_dict = find_overlap(outdir,'mobile_element',"integron")
        resistance_dict = find_overlap(outdir,'gene',"resistance")
        transposon_dict = find_overlap(outdir,'mobile_element',"transposon")
        categories = "["
        data = "["
        index = -1
        if resistance_dict != {}:
            #red
            index = index + 1
            categories = categories + "'Resistance',"
            for data_key,info in resistance_dict.items():
                name = data_key.split('|')[0]
                start = int(data_key.split('|')[1])
                end = data_key.split('|')[2]
                strand = data_key.split('|')[3]
                if len(info) == 1:
                    strand = info[0].split('|')[3]
                info_list = info[0].split('|')
                if strand == "-1":
                    start,end = end,start
                    
                    infos = info_list[0] + '|' + info_list[2] + '|' + info_list[1] + '|' + info_list[3]
                else:
                    infos = info_list[0] + '|' + info_list[1] + '|' + info_list[2] + '|' + info_list[3]
                data = data + "{name:\"%s\""%name + ",value:[%s,%s,%s,%s,\"%s\"],itemStyle:{color:'red'},label:{color:'black',rotate:45,show: true,position:'top',distance:15,textStyle:{fontSize:10},formatter:function(params){return params.name;}}},"%(index,start,end,strand,infos)
        if IS_dict != {}:
            #blue
            index_num,ret_data,categorie = config_for_each_dict(index,'blue',IS_dict,'IS')
            data = data + ret_data
            categories = categories + categorie
            index = index_num
        if integron_dict != {}:
            #yellow
            index_num,ret_data,categorie = config_for_each_dict(index,'yellow',integron_dict,'Integron')
            data = data + ret_data
            categories = categories + categorie
            index = index_num
        if transposon_dict != {}:
            #green
            index_num,ret_data,categorie = config_for_each_dict(index,'green',transposon_dict,'Transposon')
            data = data + ret_data
            categories = categories + categorie
            index = index_num
        categories = categories+"]"
        data = data + "]"
        if data and categories:
            generate_html(data,categories,endpos,outdir)


            
def exchange(start,end,strand):
    if strand == -1:
        start,end = end,start
    else:
        start,end = start,end
    return str(start),str(end)
    
def get_direction(a_dict):
    strand_list = []
    for data_key,info in a_dict.items():
        for data in info:
            each_strand = data.split('|')[-1]
            if each_strand not in strand_list:
                strand_list.append(each_strand)
    if '1' in strand_list and '-1' in strand_list:
        flag = 0
    elif '1' in strand_list and '-1' not in strand_list:
        flag = 1
    elif '1' not in strand_list and '-1' in strand_list:
        flag = -1
    return flag
    
def get_position(info):
    min = 999999999999999
    max = 0
    for data in info:
        start = int(data.split('|')[1])
        end = int(data.split('|')[2])
        if start <= min:
            min = start
        if end >= max:
            max = end
    return min,max

def config_for_each_dict(index,color,a_dict,type):
    data = ""
    categorie = ""
    flag = get_direction(a_dict)
    if flag == 0:
        index = index + 2
        categorie = categorie + "'%s (+)','%s (-)',"%(type,type)
        for data_key,info in a_dict.items():
            if (len(info)) > 1:
                start,end = get_position(info)
                name = data_key.split('|')[0]
                info_plus = []
                info_minus = []
                for datas in info:
                    each_strand1 = datas.split('|')[-1]
                    if each_strand1 == "1":
                        info_plus.append(datas)
                    else:
                        new_data = datas.split('|')[0] + '|' + datas.split('|')[2] + '|' +datas.split('|')[1] + '|-1'
                        #print(new_data)
                        info_minus.append(new_data)
                info_plus_str = ','.join(info_plus)
                info_minus_str = ','.join(info_minus)
                if info_plus_str:
                    if len(info_plus) > 1:
                        data = data + "{name:\"%s\""%name + ",value:[%s,%s,%s,%s,\"%s\"],itemStyle:{color:'%s'}},"%(index-1,str(start),str(end),'1',info_plus_str,color)
                    else:
                        data = data + "{name:\"%s\""%name + ",value:[%s,%s,%s,%s,\"%s\"],itemStyle:{color:'%s'},label:{color:'black',rotate:45,show: true,position:'top',distance:15,textStyle:{fontSize:10},formatter:function(params){return params.name;}}},"%(index-1,str(start),str(end),'1',info_plus_str,color)
                if info_minus_str:
                    if len(info_minus) > 1:
                        data = data + "{name:\"%s\""%name + ",value:[%s,%s,%s,%s,\"%s\"],itemStyle:{color:'%s'}},"%(index,str(end),str(start),'-1',info_minus_str,color)
                    else:
                        #print(info_minus_str)
                        data = data + "{name:\"%s\""%name + ",value:[%s,%s,%s,%s,\"%s\"],itemStyle:{color:'%s'},label:{color:'black',rotate:45,show: true,position:'top',distance:15,textStyle:{fontSize:10},formatter:function(params){return params.name;}}},"%(index,str(end),str(start),'-1',info_minus_str,color)
            else:
                name = info[0].split('|')[0]
                start = info[0].split('|')[1]
                end = info[0].split('|')[2]
                #infos = ','.join(info)
                each_strand = info[0].split('|')[-1]
                if each_strand == '1':
                    infos = info[0]
                    data = data + "{name:\"%s\""%name + ",value:[%s,%s,%s,%s,\"%s\"],itemStyle:{color:'%s'},label:{color:'black',rotate:45,show: true,position:'top',distance:15,textStyle:{fontSize:10},formatter:function(params){return params.name;}}},"%(index-1,start,end,each_strand,infos,color)
                else:
                    start,end = end,start
                    infos = name + '|' + start + '|' + end + '|-1'
                    data = data + "{name:\"%s\""%name + ",value:[%s,%s,%s,%s,\"%s\"],itemStyle:{color:'%s'},label:{color:'black',rotate:45,show: true,position:'top',distance:15,textStyle:{fontSize:10},formatter:function(params){return params.name;}}},"%(index,start,end,each_strand,infos,color)
    else:
        index = index + 1
        if flag == 1:
            categorie = categorie + "'%s (+)',"%type
        elif flag == -1:
            categorie = categorie + "'%s (-)',"%type
        for data_key,info in a_dict.items():
            name = data_key.split('|')[0]
            start = int(data_key.split('|')[1])
            end = data_key.split('|')[2]
            strand = data_key.split('|')[3]
            if len(info) == 1:
                strand = info[0].split('|')[3]
            if strand == "-1":
                start,end = end,start
            if (len(info)) > 1:
                if strand == "1":
                    infos = ','.join(info)
                    data = data + "{name:\"%s\""%data_key + ",value:[%s,%s,%s,%s,\"%s\"],itemStyle:{color:'%s'}},"%(index,start,end,strand,infos,color)
                else:
                    temp_info = []
                    for li in info:
                        new_li = li.split('|')[0] + '|' + li.split('|')[2] + '|' + li.split('|')[1] + '|-1'
                        temp_info.append(new_li)
                    infos = ','.join(temp_info)
                    data = data + "{name:\"%s\""%data_key + ",value:[%s,%s,%s,%s,\"%s\"],itemStyle:{color:'%s'}},"%(index,start,end,strand,infos,color)
            else:
                if strand == "1":
                    infos = info[0]
                    data = data + "{name:\"%s\""%name + ",value:[%s,%s,%s,%s,\"%s\"],itemStyle:{color:'%s'},label:{color:'black',rotate:45,show: true,position:'top',distance:15,textStyle:{fontSize:10},formatter:function(params){return params.name;}}},"%(index,start,end,strand,infos,color)
                else:
                    infos = info[0].split('|')[0] + '|' + info[0].split('|')[2] + '|' + info[0].split('|')[1] + '|-1'
                    data = data + "{name:\"%s\""%name + ",value:[%s,%s,%s,%s,\"%s\"],itemStyle:{color:'%s'},label:{color:'black',rotate:45,show: true,position:'top',distance:15,textStyle:{fontSize:10},formatter:function(params){return params.name;}}},"%(index,start,end,strand,infos,color)
    return index,data,categorie

def generate_html(data,categories,endpos,outdir):
    html_file = os.path.join(outdir,"annotaion.html")
    with open(html_file,'w') as w:
        w.write('')
    with open(html_file,'a') as w:
        w.write('''\
<!DOCTYPE html>
<html lang="en">
 <head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <meta name="Author" content="">
  <meta name="Keywords" content="">
  <meta name="Description" content="">
  <title>Graph</title>
<script type="text/javascript" src="https://cdn.jsdelivr.net/npm/echarts@4.8.0/dist/echarts.js"></script>
 </head>
 <body>
 
<div align="right">
<button id="screenshot" style="position:relative;right:2%"><a id="download">Download</a></button>
<div align="center">

<div align="center">
<div id="main" style="width:96%;height:700px" ></div>
</div>


<script type="text/javascript">

const a = window.document.querySelector('#download')
  a.addEventListener('click',()=>{
    const content = document.querySelector('svg').outerHTML
    const blob= new Blob([content], {type: 'xml/svg'})
    a.href = URL.createObjectURL(blob)
    a.download = 'Annotation.svg'
  })
var myChart = echarts.init(document.getElementById('main'),null,{renderer:'svg'});''')
        w.write("var data = %s;\nvar startpos = 0;\nvar categories = %s;\nvar endpos = %s;\n"%(data,categories,endpos))

        w.write('''
var MyShape = echarts.graphic.extendShape({
    shape: {
        x: 0,
        y: 0,
        width: 0,
        height: 0,
    },
    buildPath: function (ctx, shape) {
        if(shape.width > 0){
            ctx.moveTo(shape.x, shape.y);
            ctx.lineTo(shape.x + shape.width*15/20, shape.y);
            ctx.lineTo(shape.x + shape.width*15/20, shape.y-shape.height*1/5);
            ctx.lineTo(shape.x + shape.width, shape.y+shape.height*3/10);
            ctx.lineTo(shape.x + shape.width*15/20, shape.y+shape.height*4/5);
            ctx.lineTo(shape.x + shape.width*15/20, shape.y+shape.height*3/5);
            ctx.lineTo(shape.x , shape.y+shape.height*3/5);
            ctx.closePath();
        }
        else{
            abs_width = 0 - shape.width
            ctx.moveTo(shape.x, shape.y);
            ctx.lineTo(shape.x - abs_width*15/20, shape.y);
            ctx.lineTo(shape.x - abs_width*15/20, shape.y-shape.height*1/5);
            ctx.lineTo(shape.x - abs_width, shape.y+shape.height*3/10);
            ctx.lineTo(shape.x - abs_width*15/20, shape.y+shape.height*4/5);
            ctx.lineTo(shape.x - abs_width*15/20, shape.y+shape.height*3/5);
            ctx.lineTo(shape.x , shape.y+shape.height*3/5);
            ctx.closePath();
            }
    }
});
echarts.graphic.registerShape('myCustomShape', MyShape);

function renderItem(params, api) {
    var categoryIndex = api.value(0);
    var start = api.coord([api.value(1), categoryIndex]);
    var end = api.coord([api.value(2), categoryIndex]);
    var height = api.size([0, 1])[1] * 0.3;
    return {
        type: 'myCustomShape',
        shape: {
            x: start[0],
            y: start[1] - height / 2,
            width: end[0] - start[0],
            height: height,
        },
        style: api.style()
    };
}

option = {
    toolbox: {
        show: true,
        feature: {
            dataZoom: {
                yAxisIndex: 'none'
            },
            restore: {},
            //saveAsImage:{type:'svg'},
        }
    },
    tooltip: {
        textStyle:{
　　      align:'left'
　　　　},
        formatter: function (params) {
            var data_list = params.value[4].split(',');
            var ret_data = "";
            for (var i=0;i<data_list.length;i++){
                var name = data_list[i].split('|')[0];
                var start = data_list[i].split('|')[1];
                var end = data_list[i].split('|')[2];
                var strand = data_list[i].split('|')[3];
                if (strand == "1"){
                    ret_data = ret_data + name + ': ' + start + '-' + end + 'bp' + ' ' + '+ ; ';
                } else {
                    ret_data = ret_data + name + ': ' + end + '-' + start + 'bp' + ' ' + '- ; ';
                }
                if ((i+1)%2==0){
                    ret_data = ret_data + '<br>';
                }
            }
            return  ret_data;
        }
    },
    title: {
        text: 'Annotation',
        left: 'center',
        backgroundColor: 'transparent'
    },
    dataZoom: [{
        type: 'slider',
        filterMode: 'weakFilter',
        showDataShadow: false,
        top: 700,
        height: 20,
        borderColor: 'transparent',
        backgroundColor: '#e2e2e2',
        handleIcon: 'M10.7,11.9H9.3c-4.9,0.3-8.8,4.4-8.8,9.4c0,5,3.9,9.1,8.8,9.4h1.3c4.9-0.3,8.8-4.4,8.8-9.4C19.5,16.3,15.6,12.2,10.7,11.9z M13.3,24.4H6.7v-1.2h6.6z M13.3,22H6.7v-1.2h6.6z M13.3,19.6H6.7v-1.2h6.6z',
        handleSize: 20,
        handleStyle: {
            shadowBlur: 6,
            shadowOffsetX: 1,
            shadowOffsetY: 2,
            shadowColor: '#aaa'
        },
        labelFormatter: '',
        start: 0,
        end: 100,
    }, {
        type: 'inside',
        filterMode: 'empty',
    }],
    grid: {
        height: 600,
    },
    xAxis: {
        min: startpos,
        max: endpos,
        scale: true,
        axisLabel: {
            formatter: function (val) {
                return Math.max(0, val - startpos) + ' bp';
            }
        }
    },
    yAxis: {
        data: categories
    },
    series: [{
        type: 'custom',
        renderItem: renderItem,
        data: data,
        encode: {
            x: [1, 2],
            y: 0
        },
        
    }]
};
myChart.setOption(option);
</script>

 </body>
 
</html>''')

