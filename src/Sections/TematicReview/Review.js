import React, { Component } from 'react';
import { useState, useEffect, createRef } from 'react';

import { Navigate } from 'react-router-dom';

import { AgGridReact } from 'ag-grid-react';
import 'ag-grid-enterprise';
import 'ag-grid-community/styles/ag-grid.css';
import 'ag-grid-community/styles/ag-theme-alpine.css';

import Plot from 'react-plotly.js';

import Slider from 'react-input-slider';

import { variables } from '../Variables.js';


var topicFilterParams = {
  comparator: (TopicParam, cellValue) => {
    if (TopicParam === cellValue) {
      return 0;
    }
    if (cellValue < TopicParam) {
      return -1;
    }
    if (cellValue > TopicParam) {
      return -1;
    }
    return 0;
  },
};

export class Review extends Component {

    constructor(props) {
        super(props);

        this.gridAnaliseRef = createRef();
        this.state = {
            updateOr: false,
            loading: false,

            // Analise table
            analise_articles: [],
            analise_info: [
                {field: 'uid'},
                {field: 'titl', filter: 'agTextColumnFilter'},
                {field: 'pdat', filter: 'agTextColumnFilter'},
                {field: 'auth', filter: 'agTextColumnFilter'},
                {field: 'jour', filter: 'agTextColumnFilter'},
                {field: 'pt', filter: 'agTextColumnFilter'},
                {field: 'topic', filter: 'agNumberColumnFilter', sortable: true, filterParams: topicFilterParams},
                {field: 'prop', filter: 'agNumberColumnFilter'},
            ],
            DetailArticle: null,

            // Analise clust graph
            clust_graph: null,
            heapmap: null,
            heirarchy: null,

            // Filter topic
            current_topic: -2,
            topics: new Set(),

            //Summirise
            summarise: null,
        }
    }

    componentDidMount() {
        this.getAnalise();
        console.log('start review');
    }

    componenDidUpdate(prevState) {
      // Популярный пример (не забудьте сравнить пропсы):
      if (this.state.current_topic !== prevState.current_topic) {
        this.externalFilterChanged();
      }
    }

    onSelectionChanged = () => {
        const selectedRows = this.gridAnaliseRef.current.api.getSelectedRows();
        this.setState({DetailArticle: (selectedRows.length === 1 ? selectedRows[0] : null)})
    }

    getAnalise = (url, interval = 1000) => {
        fetch(variables.API_URL + "/api/analise/",
          {
            headers: {
                'Content-Type': 'application/json;charset=utf-8',
                'Authorization': `Token ${variables.token}`,
            },
          }
        )
        .then((res) => {
                if (res.status == 202) {
                    this.setState({loading: true})
                    setTimeout(() => {
                      return this.getAnalise(url, interval)
                    }, interval);
                }
                if (res.status == 200) {
                    return res.json()
                } else {
                    throw Error(res.statusText)
                }
            })
        .then((data) => {
          delete data.graph.layout.width;
          delete data.heapmap.layout.width;
          delete data.heirarchy.layout.width;
          var topics = new Set()
          for (let record of data.data) {
            topics.add(record.topic)
          }
          this.setState({
            analise_articles: data.data,
            DetailArticle: data.data[0],
            clust_graph: data.graph,
            heapmap: data.heapmap,
            heirarchy: data.heirarchy,
            loading: false,
            topics: [...topics],
          });
        })
        .catch((err) => {
          console.log(err);
          this.setState({data: [], dataInfo: [], DetailArticle: null, loading: false});
        });
    }

    externalFilterChanged = (newValue) => {
        this.setState({current_topic: newValue.x, summarise: null});
        this.gridAnaliseRef.current.api.onFilterChanged();
    }

    isExternalFilterPresent = () => {
        // if ageType is not everyone, then we are filtering
        return this.state.current_topic !== -2;
    }

    doesExternalFilterPass = (node) => {
      if (node.data) {
        if (node.data.topic === this.state.current_topic) {
            return true
        }
      }
      return false;
    }

    getGraphData = () => {
        var current_topic = this.state.current_topic.toString();
        console.log(current_topic)
        if (current_topic === '-2') {
            return this.state.clust_graph.data
        }

        var data = []
        for (let topic of this.state.clust_graph.data) {
            var topic_id = topic.name.split('_')[0];
            if (current_topic === topic_id) {
                data.push(topic);
            }
        }
        return data
    }

    getSummarise = (task_id, interval = 1000) => {
        fetch(variables.API_URL + `/api/summarise?task_id=${task_id}`,
          {
            headers: {
                'Content-Type': 'application/json;charset=utf-8',
                'Authorization': `Token ${variables.token}`,
            },
          }
        )
        .then((res) => {
                if (res.status == 202) {
                    this.setState({loading: true})
                    setTimeout(() => {
                      return this.getSummarise(task_id, interval)
                    }, interval);
                }
                if (res.status == 200) {
                    return res.json()
                } else {
                    throw Error(res.statusText)
                }
            })
        .then((data) => {
          this.setState({
            summarise: data.data
          });
        })
        .catch((err) => {
          console.log(err);
          this.setState({summarise: null});
        });
    }

    createSummariseQuery() {
        var data = []
        for (let article of this.state.analise_articles) {
            if (this.state.current_topic === article.topic) {
                data.push(article);
            }
        }
        fetch(variables.API_URL + '/api/summarise', {
            method: 'POST',
            headers: {
                'Accept': 'application/json',
                'Content-Type': 'application/json;charset=utf-8',
                'Authorization': `Token ${variables.token}`,
            },
            body: JSON.stringify({
                articles: data
            })
        })
            .then((res) => {
                if (res.status == 200) { return res.json() }
                else { throw Error(res.statusText) }
            })
            .then((result) => {
                var task_id = result.data;
                this.getSummarise(task_id);
            })
            .catch((error) => {
                alert('Ошибка')
            })
    }

    render() {
        const {
            token,
            analise_articles,
            analise_info,
            DetailArticle,

            clust_graph,
            heapmap,
            heirarchy,

            current_topic,
            topics,

            summarise
        } = this.state;

        return (
            <>
                            <section class="col shadow p-4" style={{backgroundColor: "#fff"}}>
                                    <div>
                                        <p>Topic {current_topic === -2? "Выбраны все": `№ ${current_topic}`}</p>
                                        <Slider
                                            axis="x"
                                            x={current_topic}
                                            xmax={topics.length - 2}
                                            xmin={-2}
                                            onChange={this.externalFilterChanged}
                                        />
                                    </div>
                                    <div className="ag-theme-alpine" style={{height: 700}}>
                                        <AgGridReact
                                            ref={this.gridAnaliseRef}
                                            rowData={analise_articles}
                                            columnDefs={analise_info}
                                            pagination={true}
                                            rowSelection={'single'}
                                            onSelectionChanged={this.onSelectionChanged}
                                            animateRows={true}
                                            isExternalFilterPresent={this.isExternalFilterPresent}
                                            doesExternalFilterPass={this.doesExternalFilterPass}
                                        >
                                        </AgGridReact>
                                    </div>
                                    <div>
                                        {clust_graph?
                                        <Plot
                                            data={this.getGraphData()}
                                            layout={clust_graph.layout}
                                          />
                                        :null
                                        }
                                    </div>
                                    <div>
                                        {heapmap?
                                        <Plot
                                            data={heapmap.data}
                                            layout={heapmap.layout}
                                          />
                                        :null
                                        }
                                    </div>
                                    <div>
                                        {heirarchy?
                                        <Plot
                                            data={heirarchy.data}
                                            layout={heirarchy.layout}
                                          />
                                        :null
                                        }
                                    </div>
                                    <div>
                                        {summarise?
                                        <>
                                            <p>Summarise</p>
                                            <p>{summarise}</p>
                                        </>
                                        :<input className="btn btn-primary" type="submit" value="Суммаризовать" onClick={() => this.createSummariseQuery()}/>}
                                    </div>
                            </section>

                            <aside id="sidebar2" class="col-md-4 bg-light collapse show width mb-5 shadow">
                              <h1 class="h2 pt-3 pb-2 mb-3 border-bottom">Подробности</h1>
                              <nav class="small" id="toc">
                                {DetailArticle?
                                    <div class="card mb-3">
                                        <div class="card-body">
                                          <a href= { DetailArticle.url } class="card-title link-primary text-decoration-none h5"> { DetailArticle.titl } </a>
                                          <p class="card-text">---------------------------------- </p>
                                          <p class="card-text">Авторы :  { DetailArticle.auth } </p>
                                          <p class="card-text">---------------------------------- </p>
                                          <p class="card-text">Аннотация :  </p>
                                          <p class="card-text"> { DetailArticle.tiab } </p>
                                          <p class="card-text">---------------------------------- </p>
                                          <p class="card-text"><small class="text-success">Дата публикации : { DetailArticle.pdat } </small></p>
                                          <p class="card-text"><small class="text-success">Издание : { DetailArticle.jour }</small></p>
                                          <p class="card-text"><small class="text-success">Вид публикации : { DetailArticle.pt }</small></p>
                                          <p class="card-text"><small class="text-success">Страна : { DetailArticle.pl } </small></p>
                                          <p class="card-text"><small class="text-success">{ DetailArticle.mesh } </small></p>
                                        </div>
                                      </div>
                                :null}
                              </nav>
                            </aside>
            </>
        )
    }
}