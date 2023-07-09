import React, { Component } from 'react';
import { useState, useEffect, createRef } from 'react';

import { Navigate } from 'react-router-dom';

import Graph from "react-graph-vis";
import { v4 as uuidv4 } from 'uuid'
//import "./network.css";

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

export class TematicReview extends Component {

    constructor(props) {
        super(props);

        this.gridRef = createRef();
        this.gridAnaliseRef = createRef();
        this.state = {
            updateOr: false,
            loading: false,
            useAll: true,
            token: variables.token,

        // Search
            articles: [],
            DetailArticle: null,
            articlesInfo: [
                {
                    field: 'uid',

                },
                {field: 'titl', filter: 'agTextColumnFilter'},
                {field: 'pdat', filter: 'agTextColumnFilter'},
                {field: 'auth', filter: 'agTextColumnFilter'},
                {field: 'jour', filter: 'agTextColumnFilter'},
                {field: 'pt', filter: 'agTextColumnFilter'},
                {field: 'topic', filter: 'agNumberColumnFilter', sortable: true, filterParams: topicFilterParams},
                {field: 'prop', filter: 'agNumberColumnFilter'},
            ],
            translation_stack: null,
            full_query: null,
            short_query: null,
            task: null,
            message: null,
            count: 0,
            analiseRows: [],

            // Filters
            queryText: 'covid-19',
            queryStartDate: '2022-01-01',
            queryEndDate: new Date().toISOString().split('T')[0],
            queryTypes: new Set(),
            queryOlds: new Set(),
            queryGenders: new Set(),

        // Analise
            // Analise table
            analise_articles: [],
            analise_info: [
                {field: 'uid'},
                {field: 'titl', filter: 'agTextColumnFilter'},
                {field: 'pdat', filter: 'agTextColumnFilter'},
                {field: 'auth', filter: 'agTextColumnFilter'},
                {field: 'jour', filter: 'agTextColumnFilter'},
                {field: 'pt', filter: 'agTextColumnFilter'},
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

        //Other
            graph: null,
        }
    }


    componentDidMount() {
        this.getArticles();
        console.log('start');
    }

    // Search

    createQuery() {
        let query = "/api"
        if (this.state.queryText) {
            query = `${query}?search_field=${this.state.queryText}`
        }
        if (this.state.queryStartDate) {
            query = `${query}&dateStart=${this.state.queryStartDate}`
        }
        if (this.state.queryEndDate) {
            query = `${query}&dateStop=${this.state.queryEndDate}`
        }
        if (this.state.queryGenders.size != 0) {
            for (let el of [...this.state.queryGenders]) {
                query = `${query}&Gender=${el}`
            }
        }
        if (this.state.queryTypes.size != 0) {
            for (let el of [...this.state.queryTypes]) {
                query = `${query}&Type=${el}`
            }
        }
        if (this.state.queryOlds.size != 0) {
            for (let el of [...this.state.queryOlds]) {
                query = `${query}&Old=${el}`
            }
        }

        if (query === "/api") {
            throw "Not query text"
        }

        return query
    }

//    getArticles = (url, interval = 1000) => {
//      fetch(variables.API_URL + '/api/articles/',
//            {
//                headers: {
//                    'Content-Type': 'application/json;charset=utf-8',
//                    'Authorization': `Token ${variables.token}`,
//                },
//            }
//          )
//          .then(response => {
//                console.log(response.status);
//                if (response.status == 202) {
//                    this.setState({loading: true})
//                    setTimeout(() => {
//                      return this.getArticles(url, interval)
//                    }, interval);
//                }
//                if (response.status == 200) {
//                    return response.json()
//                } else {
//                    throw Error(response.statusText)
//                }
//          })
//          .then(data => {
//            this.setState({
//                articles: data.data, DetailArticle: data.data[0], message: data.message, loading: false
//            });
//          })
//          .catch(error => {
//            console.log(error);
//            this.setState({ articles: [], articlesInfo: [], loading: false });
//          })
//    }

    getArticles = (url, interval = 1000) => {
      fetch(variables.API_URL + '/api/articles/',
            {
                headers: {
                    'Content-Type': 'application/json;charset=utf-8',
                    'Authorization': `Token ${variables.token}`,
                },
            }
          )
          .then(response => {
                console.log(response.status);
                if (response.ok) {
                    return response.json()
                } else {
                    throw Error(response.statusText)
                }
          })
          .then(data => {
            if (!data.data) {
                this.setState({loading: true, message: data.message});
                console.log(data.message);
                setTimeout(() => {
                  return this.getArticles(url, interval)
                }, interval);
            } else {
                this.setState({
                    articles: data.data, DetailArticle: data.data[0], message: data.message, loading: false
                });
            }
          })
          .catch(error => {
            console.log(error);
            this.setState({ articles: [], articlesInfo: [], loading: false, message: 'Something gone wrong' });
          })
    }

    createTask() {
        // Отправляем запрос на сервер для получения статей
        let query_url = ''
        try {
            query_url = variables.API_URL + this.createQuery()
        } catch(e) {
            alert('Не введен текст для поиска');
            return
        }
        fetch(query_url,
            {
                headers: {
                    'Content-Type': 'application/json;charset=utf-8',
                    'Authorization': `Token ${this.state.token}`,
                },
            })
            .then(response => response.json())
            .then(data => {
                this.setState({
                    full_query: data.task.full_query,
                    translation_stack: data.task.translation_stack,
                    short_query: data.task.query,
                    task: data.task,
                    count: data.task.count,
                    articles: [],
                    articlesInfo: [],
                    loading: true,

                });
                alert("You query in queue? please wait to get result");
                this.getArticles()
            })
            .catch(error => {
                console.log(error);
                this.setState({ task: null });
            }
        )
    }

    componentDidMount() {
        this.getArticles();
        this.getAnalise();
        console.log('start search');
    }

    onSelectionAnalise = () => {
        const selectedRows = this.gridRef.current.api.getSelectedRows();
        this.setState({DetailArticle: (selectedRows.length === 1 ? selectedRows[0] : null)})
    }

    changeQueryText = (e) => {
        this.setState({ queryText: e.target.value });
    }

    changeQueryStartDate = (e) => {
        this.setState({ queryStartDate: e.target.value });
    }

    changeQueryEndDate = (e) => {
        this.setState({ queryEndDate: e.target.value });
    }

    changeQueryGenders(gender) {
        if (this.state.queryGenders.has(gender)) {
            this.state.queryGenders.delete(gender)
        } else {
            this.state.queryGenders.add(gender)
        }
        this.setState({updateOr: !this.state.updateOr})

    }

    changeQueryTypes(type) {
        if (this.state.queryTypes.has(type)) {
            this.state.queryTypes.delete(type)
        } else {
            this.state.queryTypes.add(type)
        }
        this.setState({updateOr: !this.state.updateOr})
    }

    changeQueryOlds(old) {
        if (this.state.queryOlds.has(old)) {
            this.state.queryOlds.delete(old)
        } else {
            this.state.queryOlds.add(old)
        }
        this.setState({updateOr: !this.state.updateOr})
    }

    changeQueryOldsMany(olds) {
        if (!this.state.useAll) {
            for (let old of olds) {
                this.state.queryOlds.delete(old)
            }
        } else {
            for (let old of olds) {
                this.state.queryOlds.add(old)
            }
        }
        this.setState({useAll: !this.state.useAll})
    }

    RoundPersent(number) {
        return number.toFixed(2);
    }

    startAnalise() {
        let analise_data = [];
        this.gridRef.current.api.forEachNodeAfterFilter((rowNode) => analise_data.push(rowNode.data.uid));
        fetch(variables.API_URL + '/api/analise/', {
            method: 'POST',
            headers: {
                'Accept': 'application/json',
                'Content-Type': 'application/json;charset=utf-8',
                'Authorization': `Token ${variables.token}`,
            },
            body: JSON.stringify({
                articles: analise_data
            })
        })
            .then((res) => {
                if (res.status == 200) { return res.json() }
                else { throw Error(res.statusText) }
            })
            .then((result) => {
                this.getAnalise();
            })
            .catch((error) => {
                alert('Ошибка')
            })
    }

    // Analise

    onSelectionChanged = (gridApi) => {
        const selectedRows = this.gridAnaliseRef.current.api.getSelectedRows();
        this.setState({DetailArticle: (selectedRows.length === 1 ? selectedRows[0] : null)})
    }

//    getAnalise = (url, interval = 1000) => {
//        fetch(variables.API_URL + "/api/analise/",
//          {
//            headers: {
//                'Content-Type': 'application/json;charset=utf-8',
//                'Authorization': `Token ${variables.token}`,
//            },
//          }
//        )
//        .then((res) => {
//                if (res.status == 202) {
//                    this.setState({loading: true})
//                    setTimeout(() => {
//                      return this.getAnalise(url, interval)
//                    }, interval);
//                }
//                if (res.status == 200) {
//                    return res.json()
//                } else {
//                    throw Error(res.statusText)
//                }
//            })
//        .then((data) => {
//          delete data.graph.layout.width;
//          delete data.heapmap.layout.width;
//          delete data.heirarchy.layout.width;
//          var topics = new Set()
//          for (let record of data.data) {
//            topics.add(record.topic)
//          }
//          this.setState({
//            analise_articles: data.data,
//            DetailArticle: data.data[0],
//            clust_graph: data.graph,
//            heapmap: data.heapmap,
//            heirarchy: data.heirarchy,
//            loading: false,
//            topics: [...topics],
//          });
//        })
//        .catch((err) => {
//          console.log(err);
//          this.setState({data: [], dataInfo: [], DetailArticle: null, loading: false});
//        });
//    }

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
                if (res.ok) {
                    return res.json()
                } else {
                    throw Error(res.statusText)
                }
            })
        .then((data) => {
          if (!data.data) {
            this.setState({loading: true, message: data.message});
                console.log(data.message);
                setTimeout(() => {
                  return this.getAnalise(url, interval)
                }, interval);
          } else {
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
                message: data.message,
                topics: [...topics],
              });
          }
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
            count,
            articles,
            DetailArticle,
            articlesInfo,
            translation_stack,
            full_query,
            short_query,
            message,
            loading,

            queryText,
            queryEndDate,
            queryStartDate,

            analise_articles,
            analise_info,
            clust_graph,
            heapmap,
            heirarchy,
            current_topic,
            topics,
            summarise,

            graph
        } = this.state;

        if (!token){
            return <Navigate push to="/login" />
        } else {
            return (
            <>
                <div className="container-fluid">
                    <div className="row">
                      <header className="navbar navbar-dark sticky-top bg-primary flex-md-nowrap p-0 shadow-sm">
                      <div className="col-md-2">
                        <div className="row g-0">
                          <div className="col-md-2">
                            <button className="mt-1 navbar-toggler collapsed" type="button" data-bs-toggle="collapse" data-bs-target="#sidebar" aria-controls="sidebar" aria-expanded="false" aria-label="Toggle navigation">
                              <span className="icon-bar top-bar"></span>
                              <span className="icon-bar middle-bar"></span>
                              <span className="icon-bar bottom-bar"></span>
                            </button>
                          </div>
                          <div className="col-md-10">
                            <a className="navbar-brand" href="#">SECHENOV AI-DATAMED</a>
                          </div>
                        </div>
                      </div>
                      <div className="col-md-8">
                        <input
                            className="form-control w-100"
                            id="search"
                            type="text"
                            name="search_field"
                            placeholder="Поисковый запрос"
                            value={queryText}
                            onChange={this.changeQueryText}
                            aria-label="Search" />
                      </div>
                      <div className="col-md-2">
                        <div className="row g-0">
                          <div className="col-md-10">
                            <ul className="navbar-nav px-3">
                              <li className="nav-item text-nowrap">
                                <input className="btn btn-primary" type="submit" value="Найти" onClick={() => this.createTask()}/>
                              </li>
                            </ul>
                          </div>
                          <div className="col-md-2">
                            <button className="mt-2 navbar-toggler collapsed" type="button" data-bs-toggle="collapse" data-bs-target="#sidebar2" aria-controls="sidebar2" aria-expanded="false" aria-label="Toggle navigation">
                              <span className="icon-bar top-bar"></span>
                              <span className="icon-bar middle-bar"></span>
                              <span className="icon-bar bottom-bar"></span>
                            </button>
                          </div>
                        </div>
                      </div>
                      </header>
                    </div>
                </div>
                <main>
                    <div>
                        <div className="container-fluid">
                            <div className="row align-items-stretch b-height">

                            <aside id="sidebar" className="col-md-2 bg-light collapse show width mb-5 shadow-sm g-0">
                                <div className="accordion accordion-flush" id="accordionFlushExample">
                                  <div className="accordion-item">
                                    <h2 className="accordion-header" id="flush-headingOne">
                                      <button className="accordion-button collapsed" type="button" data-bs-toggle="collapse" data-bs-target="#flush-collapseOne" aria-expanded="false" aria-controls="flush-collapseOne">
                                        Дата публикации
                                      </button>
                                    </h2>
                                    <div id="flush-collapseOne" className="collapse show multi-collapse" aria-labelledby="flush-headingOne" data-bs-target="#accordionFlushExample">
                                      <div className="accordion-body">
                                        <div className="mb-3">
                                          <label for="localdate">От : </label>
                                          <input
                                            type="date"
                                            id="d1"
                                            name="dateStart"
                                            value={queryStartDate}
                                            onChange={this.changeQueryStartDate}
                                          />
                                        </div>
                                        <div>
                                          <label for="localdate">До : </label>
                                          <input
                                            type="date"
                                            id="d2"
                                            name="dateStop"
                                            value={queryEndDate}
                                            onChange={this.changeQueryEndDate}
                                          />
                                        </div>
                                      </div>
                                    </div>
                                  </div>
                                  <div className="accordion-item">
                                    <h2 class="accordion-header" id="flush-headingThree">
                                      <button class="accordion-button collapsed" data-target='#flush-collapseThree' type="button" data-bs-toggle="collapse" data-bs-target="#flush-collapseFour" aria-expanded="false" aria-controls="flush-collapseThree">
                                        Тип статьи
                                      </button>
                                    </h2>
                                    <div id="flush-collapseFour" className="collapse show multi-collapse" aria-labelledby="flush-headingFour" data-bs-target="#accordionFlushExample">
                                      <div className="accordion-body">
                                        <div className="form-check form-check-inline">
                                          <input
                                            className="form-check-input"
                                            type="checkbox"
                                            id="CheckboxClinicalTrial"
                                            name = "CheckBoxClinicalTrial"
                                            checked={this.state.queryTypes.has('clinical trial')}
                                            onChange={() => this.changeQueryTypes('clinical trial')}
                                          />
                                          <label className="form-check-label" for="CheckboxClinicalTrial">Clinical Trial</label>
                                        </div>
                                        <div className="form-check form-check-inline">
                                          <input
                                            className="form-check-input"
                                            type="checkbox"
                                            id="CheckboxMetaAnalysys"
                                            name = "CheckboxMetaAnalysys"
                                            checked={this.state.queryTypes.has('meta-analysis')}
                                            onChange={() => this.changeQueryTypes('meta-analysis')}
                                          />
                                          <label className="form-check-label" for="CheckboxMetaAnalysys">Meta Analysys</label>
                                        </div>
                                        <div className="form-check form-check-inline">
                                          <input
                                            className="form-check-input"
                                            type="checkbox"
                                            id="CheckboxRandomizedControlledTrial"
                                            name ="CheckboxRandomizedControlledTrial"
                                            checked={this.state.queryTypes.has('randomized controlled trial')}
                                            onChange={() => this.changeQueryTypes('randomized controlled trial')}
                                          />
                                          <label className="form-check-label" for="CheckboxRandomizedControlledTrial">Randomized Controlled Trial</label>
                                        </div>
                                        <div className="form-check form-check-inline">
                                          <input
                                            className="form-check-input"
                                            type="checkbox"
                                            id="CheckboxReview"
                                            name ="CheckboxReview"
                                            checked={this.state.queryTypes.has('review')}
                                            onChange={() => this.changeQueryTypes('review')}
                                          />
                                          <label className="form-check-label" for="CheckboxReview">Review</label>
                                        </div>
                                        <div className="form-check form-check-inline">
                                          <input
                                            className="form-check-input"
                                            type="checkbox"
                                            id="CheckboxSystematicReview"
                                            name ="CheckboxSystematicReview"
                                            checked={this.state.queryTypes.has('systematic review')}
                                            onChange={() => this.changeQueryTypes('systematic review')}
                                          />
                                          <label className="form-check-label" for="CheckboxSystematicReview">Systematic Review</label>
                                        </div>
                                        <div className="form-check form-check-inline">
                                          <input
                                            className="form-check-input"
                                            type="checkbox"
                                            id="CheckboxJournalArticle"
                                            name ="CheckboxJournalArticle"
                                            checked={this.state.queryTypes.has('journal article')}
                                            onChange={() => this.changeQueryTypes('journal article')}
                                          />
                                          <label className="form-check-label" for="CheckboxJournalArticle">Journal Article</label>
                                        </div>

                                      </div>
                                    </div>
                                  </div>
                                  <div className="accordion-item">
                                    <h2 className="accordion-header" id="flush-headingThree">
                                      <button className="accordion-button collapsed" data-target='#flush-collapseThree' type="button" data-bs-toggle="collapse" data-bs-target="#flush-collapseThree" aria-expanded="false" aria-controls="flush-collapseThree">
                                        Возраст пациента
                                      </button>
                                    </h2>
                                    <div id="flush-collapseThree" className="collapse show multi-collapse" aria-labelledby="flush-headingThree" data-bs-target="#accordionFlushExample">
                                      <div className="accordion-body">
                                        <div className="form-check form-check-inline">
                                          <input
                                            className="form-check-input"
                                            type="checkbox"
                                            id="CheckboxChild18"
                                            checked={this.state.queryOlds.has('infant[mh]') & this.state.queryOlds.has('child[mh]') & this.state.queryOlds.has('adolescent[mh]')}
                                            onChange={() => this.changeQueryOldsMany(['infant[mh]','child[mh]', 'adolescent[mh]'])}
                                          />
                                          <label className="form-check-label" for="CheckboxChild18">Child: birth-18 years</label>
                                        </div>
                                        <div className="form-check form-check-inline">
                                          <input
                                            className="form-check-input"
                                            type="checkbox"
                                            id="CheckboxNewborn"
                                            checked={this.state.queryOlds.has('infant, newborn[mh]')}
                                            onChange={() => this.changeQueryOlds('infant, newborn[mh]')}
                                          />
                                          <label className="form-check-label" for="CheckboxNewborn">Newborn: birth-1 months</label>
                                        </div>
                                        <div className="form-check form-check-inline">
                                          <input
                                            className="form-check-input"
                                            type="checkbox"
                                            id="CheckboxInfant023"
                                            checked={this.state.queryOlds.has('infant[mh]')}
                                            onChange={() => this.changeQueryOlds('infant[mh]')}
                                          />
                                          <label className="form-check-label" for="CheckboxInfant023">Infant: birth-23 months</label>
                                        </div>
                                        <div className="form-check form-check-inline">
                                          <input
                                            className="form-check-input"
                                            type="checkbox"
                                            id="CheckboxInfant123"
                                            checked={this.state.queryOlds.has('infant[mh:noexp]')}
                                            onChange={() => this.changeQueryOlds('infant[mh:noexp]')}
                                          />
                                          <label className="form-check-label" for="CheckboxInfant123">Infant: 1-23 months</label>
                                        </div>
                                        <div className="form-check form-check-inline">
                                          <input
                                            className="form-check-input"
                                            type="checkbox"
                                            id="CheckboxPreschool"
                                            checked={this.state.queryOlds.has('child, preschool[mh]')}
                                            onChange={() => this.changeQueryOlds('child, preschool[mh]')}
                                          />
                                          <label className="form-check-label" for="CheckboxPreschool">Preschool Child: 2-5 years</label>
                                        </div>
                                        <div className="form-check form-check-inline">
                                          <input
                                            className="form-check-input"
                                            type="checkbox"
                                            id="CheckboxChild612"
                                            checked={this.state.queryOlds.has('child[mh:noexp]')}
                                            onChange={() => this.changeQueryOlds('child[mh:noexp]')}
                                          />
                                          <label className="form-check-label" for="CheckboxChild612">Child: 6-12 years</label>
                                        </div>
                                        <div className="form-check form-check-inline">
                                          <input
                                            className="form-check-input"
                                            type="checkbox"
                                            id="CheckboxAdolescent"
                                            checked={this.state.queryOlds.has('adolescent[mh]')}
                                            onChange={() => this.changeQueryOlds('adolescent[mh]')}
                                          />
                                          <label className="form-check-label" for="CheckboxAdolescent">Adolescent: 13-18 years</label>
                                        </div>
                                        <div className="form-check form-check-inline">
                                          <input
                                            className="form-check-input"
                                            type="checkbox"
                                            id="CheckboxAdult19"
                                            checked={this.state.queryOlds.has('adult[mh]')}
                                            onChange={() => this.changeQueryOlds('adult[mh]')}
                                          />
                                          <label className="form-check-label" for="CheckboxAdult19">Adult: 19+ years</label>
                                        </div>
                                        <div className="form-check form-check-inline">
                                          <input
                                            className="form-check-input"
                                            type="checkbox"
                                            id="CheckboxYoungAdult"
                                            checked={this.state.queryOlds.has('young adult[mh]')}
                                            onChange={() => this.changeQueryOlds('young adult[mh]')}
                                          />
                                          <label className="form-check-label" for="CheckboxYoungAdult">Young Adult: 19-24 years</label>
                                        </div>
                                        <div className="form-check form-check-inline">
                                          <input
                                            className="form-check-input"
                                            type="checkbox"
                                            id="CheckboxAdult1944"
                                            checked={this.state.queryOlds.has('adult[mh:noexp]')}
                                            onChange={() => this.changeQueryOlds('adult[mh:noexp]')}
                                          />
                                          <label className="form-check-label" for="CheckboxAdult1944">Adult: 19-44 years</label>
                                        </div>
                                        <div className="form-check form-check-inline">
                                          <input
                                            className="form-check-input"
                                            type="checkbox"
                                            id="CheckboxMiddle45"
                                            checked={this.state.queryOlds.has('middle aged[mh]') & this.state.queryOlds.has('aged[mh]')}
                                            onChange={() => this.changeQueryOldsMany(['middle aged[mh]', 'aged[mh]'])}
                                          />
                                          <label className="form-check-label" for="CheckboxMiddle45">Middle Aged: + Aged 45+ years</label>
                                        </div>
                                        <div className="form-check form-check-inline">
                                          <input
                                            className="form-check-input"
                                            type="checkbox"
                                            id="CheckboxMiddle4565"
                                            checked={this.state.queryOlds.has('middle aged[mh]')}
                                            onChange={() => this.changeQueryOlds('middle aged[mh]')}
                                          />
                                          <label className="form-check-label" for="CheckboxMiddle4565">Middle Aged: 45-65 years</label>
                                        </div>
                                        <div className="form-check form-check-inline">
                                          <input
                                            className="form-check-input"
                                            type="checkbox"
                                            id="CheckboxAged"
                                            checked={this.state.queryOlds.has('aged[mh]')}
                                            onChange={() => this.changeQueryOlds('aged[mh]')}
                                          />
                                          <label className="form-check-label" for="CheckboxAged">Aged: 65+ years</label>
                                        </div>
                                        <div className="form-check form-check-inline">
                                          <input
                                            className="form-check-input"
                                            type="checkbox"
                                            id="Checkbox80"
                                            checked={this.state.queryOlds.has('aged, 80 and over[mh]')}
                                            onChange={() => this.changeQueryOlds('aged, 80 and over[mh]')}
                                          />
                                          <label className="form-check-label" for="Checkbox80">80 and over: 80+ years</label>
                                        </div>
                                      </div>
                                    </div>
                                  </div>
                                    <div className="accordion-item">
                                    <h2 className="accordion-header" id="flush-headingTwo">
                                      <button className="accordion-button collapsed" type="button" data-bs-toggle="collapse" data-bs-target="#flush-collapseTwo" aria-expanded="false" aria-controls="flush-collapseTwo">
                                        Пол
                                      </button>
                                    </h2>
                                    <div id="flush-collapseTwo" className="collapse show multi-collapse" aria-labelledby="flush-headingTwo" data-bs-target="#accordionFlushExample">
                                      <div className="accordion-body">
                                        <div className="form-check form-check-inline">
                                          <input
                                            type="checkbox"
                                            id="CheckboxMan"
                                            name = "CheckboxMan"
                                            className="form-check-input"
                                            checked={this.state.queryGenders.has('man')}
                                            onChange={() => this.changeQueryGenders('man')}
                                          />
                                          <label className="form-check-label" for="CheckboxMan">Мужчина</label>
                                        </div>
                                        <div className="form-check form-check-inline">
                                          <input
                                            className="form-check-input"
                                            type="checkbox"
                                            name = "CheckboxWomen"
                                            id="CheckboxWomen"
                                            checked={this.state.queryGenders.has('woman')}
                                            onChange={() => this.changeQueryGenders('woman')}
                                          />
                                          <label className="form-check-label" for="CheckboxWomen">Женщина</label>
                                        </div>
                                      </div>
                                    </div>
                                  </div>
                                </div>
                            </aside>

                            <section class="col shadow p-4" style={{backgroundColor: "#fff"}}>
                              <div class="accordion accordion-flush" id="accordion">
                                <div class="accordion-item">
                                  <h2 class="accordion-header" id="">
                                    <button class="accordion-button collapsed" type="button" data-bs-toggle="collapse" data-bs-target="#flush-collapseSeven" aria-expanded="false" aria-controls="flush-collapseSeven">
                                      По запросу найдено { count } источников. {message}
                                    </button>
                                  </h2>
                                  <div id="flush-collapseSeven" class="collapse multi-collapse" aria-labelledby="flush-headingSeven" data-bs-target="#accordionFlushExample">
                                    <div class="accordion-body">
                                      <p class="pb-2 mb-3 border-bottom"> Запрос {short_query} .</p>
                                      <p class="pb-2 mb-3 border-bottom"> Запрос автоматически расширен до следующего вида - {full_query}.</p>
                                      <p class="pb-2 mb-3 border-bottom"> Служебная информация для анализа : {translation_stack}.</p>
                                    </div>
                                  </div>
                                </div>
                              </div>
                              <div>
                              <div class="bd-example">
                                <ul class="nav nav-pills mb-3" id="myTab" role="tablist">
                                  <li class="nav-item" role="presentation">
                                    <button class="nav-link active" id="home-tab" data-bs-toggle="tab" data-bs-target="#home" type="button" role="tab" aria-controls="home" aria-selected="false">Результаты поиска</button>
                                  </li>
                                  <li class="nav-item" role="presentation">
                                    <button class="nav-link" id="profile-tab" data-bs-toggle="tab" data-bs-target="#profile" type="button" role="tab" aria-controls="profile" aria-selected="true">Тематическое описание коллекции</button>
                                  </li>
                                  <li class="nav-item" role="presentation">
                                    <button class="nav-link" id="contact-tab" data-bs-toggle="tab" data-bs-target="#contact" type="button" role="tab" aria-controls="contact" aria-selected="false" >Схема</button>
                                  </li>
                                </ul>
                                <div class="tab-content" id="myTabContent">
                                  <div class="tab-pane fade active show" id="home" role="tabpanel" aria-labelledby="home-tab">
                                    <div class="container-fluid g-0">
                                        <div className="ag-theme-alpine" style={{height: 700}}>
                                            <AgGridReact
                                                ref={this.gridRef}
                                                rowData={articles}
                                                columnDefs={articlesInfo}
                                                pagination={true}
                                                onSelectionChanged={this.onSelectionAnalise}
                                                rowSelection={'single'}
                                                sideBar={{
                                                  toolPanels: [
                                                    {
                                                      id: 'columns',
                                                      labelDefault: 'Columns',
                                                      labelKey: 'columns',
                                                      iconKey: 'columns',
                                                      toolPanel: 'agColumnsToolPanel',
                                                      minWidth: 225,
                                                      width: 225,
                                                      maxWidth: 225,
                                                    },
                                                    {
                                                      id: 'filters',
                                                      labelDefault: 'Filters',
                                                      labelKey: 'filters',
                                                      iconKey: 'filter',
                                                      toolPanel: 'agFiltersToolPanel',
                                                      minWidth: 180,
                                                      maxWidth: 400,
                                                      width: 250,
                                                    },
                                                  ],
                                                  position: 'left',
                                                  defaultToolPanel: 'filters',
                                                }}
                                            >
                                            </AgGridReact>
                                        </div>
                                        <div>
                                            <input className="btn btn-primary" type="submit" value="Начать обработку" onClick={() => this.startAnalise()}/>
                                        </div>
                                    </div>
                                  </div>
                                  <div class="tab-pane fade" id="profile" role="tabpanel" aria-labelledby="profile-tab">
                                    <div>
                                        <p>Topic {current_topic === -2? "Выбраны все": `№ ${current_topic}`}</p>
                                        <Slider
                                            axis="x"
                                            x={current_topic}
                                            xmax={topics.length - 1}
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
                                            sideBar={{
                                                  toolPanels: [
                                                    {
                                                      id: 'columns',
                                                      labelDefault: 'Columns',
                                                      labelKey: 'columns',
                                                      iconKey: 'columns',
                                                      toolPanel: 'agColumnsToolPanel',
                                                      minWidth: 225,
                                                      width: 225,
                                                      maxWidth: 225,
                                                    },
                                                    {
                                                      id: 'filters',
                                                      labelDefault: 'Filters',
                                                      labelKey: 'filters',
                                                      iconKey: 'filter',
                                                      toolPanel: 'agFiltersToolPanel',
                                                      minWidth: 180,
                                                      maxWidth: 400,
                                                      width: 250,
                                                    },
                                                  ],
                                                  position: 'left',
                                                  defaultToolPanel: 'filters',
                                                }}
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
                                  </div>
                                  <div class="tab-pane fade" id="contact" role="tabpanel" aria-labelledby="contact-tab">
                                    <div class="container-fluid g-0">
                                      <div id="mynetwork" style={{width: "100%", height: "600px"}}>
                                      {graph?
                                        <Graph
                                          graph={graph}
                                          zoomView={true}
//                                          events={this.selectEvent}
                                          key={uuidv4()}
                                          options={{
                                            layout: {
                                                hierarchical: false,
                                            },
                                            edges: {
                                              color: "#000000"
                                            }
                                          }}
//                                          getNetwork={network => {
//                                            setNetwork(network);
//                                          }}
                                        />
                                        :null}
                                      </div>
                                    </div>
                                  </div>
                                </div>
                              </div>
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

                            </div>
                        </div>
                    </div>
                </main>
            </>
            )
        }
    }
}