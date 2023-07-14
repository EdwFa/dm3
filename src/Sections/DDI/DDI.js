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
import '../ag-theme-acmecorp.css';


import { variables } from '../Variables.js';

import Slider from 'react-input-slider';


const obj_color = {
    'disease': "#fdbbbb",
    'drug': "#ECC58B",
    'gene': "#E2DB8C",
    'chemical': "#21c354",
    'species': "#A6EFDC",
    'mutation': "#B2DDEA",
    'cell_type': "#C6DEF5",
    'cell_line': "#A3B3D2",
    'DNA': "#C9B9E8",
    'RNA': "#D7DBE8",
}

function markup_text(text, annotations) {
    if (!annotations) {
        return text
    }
    let markup_text = ''
    let last_position = 0
    for (let annotation of annotations) {
        let start = annotation.span.begin
        let end = annotation.span.end
        if (!annotation.prop) { console.log('This')}
        let markup_str = `<span style=\"color: ${obj_color[annotation.obj]}\">${text.slice(start, end)}<sub>${annotation.prob? annotation.prob.toFixed(2): ""}</sub></span>`
        markup_text = `${markup_text}${text.slice(last_position, start)}${markup_str}`

        last_position = end
    }
    return markup_text
}


export class DDIReview extends Component {

    constructor(props) {
        super(props);

        this.gridRef = createRef();
        this.gridAnaliseRef = createRef();
        this.state = {
            token: variables.token,
            loading: false,

            //queries
            query_list: [],
            message: null,
            articles: [],
            articlesInfo: [
                {field: 'query_number'},
                {field: 'score', filter: 'agNumberColumnFilter', 'sortable': true},
                {field: 'text', filter: 'agTextColumnFilter'},
                {field: 'section', filter: 'agTextColumnFilter'},
            ],
            summarise: null,
            task_id: null,

            // Filters
            queryText: 'What methods are available to measure anti-mullerian hormone concentrations in young women?',
            queryStartDate: '2022-01-01',
            queryEndDate: null,
            queryScore: 0.8,
            queryTypes: new Set(),
        }
    }

    getArticles = (task_id, query_number = 0, interval = 1000) => {
      fetch(variables.API_URL + `/api/ddi_review`,
            {
                headers: {
                    'Content-Type': 'application/json;charset=utf-8',
                    'Authorization': `Token ${variables.token}`,
                },
            }
          )
          .then(response => {
                console.log(query_number)
                console.log(response.status);
                if (response.ok) {
                    return response.json()
                } else {
                    throw Error(response.statusText)
                }
          })
          .then(data => {
            if (data.data === null) {
                this.setState({loading: true, message: data.message});
                setTimeout(() => {
                  return this.getArticles(task_id, query_number, interval)
                }, interval);
            } else {
                this.setState({
                    articles: [...this.state.articles, ...data.data], DetailArticle: data.data[0], loading: false,  message: 'Выполненно'
                });
                if (query_number !== 0) {
                    this.state.query_list[query_number - 1].status = 1;
                }
            }
          })
          .catch(error => {
            console.log(error);
            this.setState({ articles: [], DetailArticle: null, loading: false,  message: 'Что-то пошло не так' });
            if (query_number !== 0) {
                    this.state.query_list[query_number - 1].status = 2;
                }
          })
    }

    createTask() {
        // Отправляем запрос на сервер для получения статей
        this.state.query_list.push({query: this.state.queryText, status: 0});
        const query_number = this.state.query_list.length;
        fetch(variables.API_URL + '/api/ddi_review', {
            method: 'POST',
            headers: {
                'Accept': 'application/json',
                'Content-Type': 'application/json;charset=utf-8',
                'Authorization': `Token ${variables.token}`,
            },
            body: JSON.stringify({
                query: this.state.queryText,
                score: this.state.queryScore,
                number_of_query: query_number
            })
            })
            .then(response => response.json())
            .then(data => {
                this.setState({
                    task_id: data.data,
                    message: 'Запрос начал обрабатываться'
                });
                alert("Ваш запрос в очереди. Пожайлуста дождитесь результата");
                this.getArticles(data.data, query_number)
            })
            .catch(error => {
                console.log(error);
                this.setState({ task: null, message: 'Что-то пошло не так' });
            }
        )
    }

    clearTask() {
        this.setState({query_list: [], articles: [], DetailArticle: null})
        alert("Таблица очищена!");
    }

    componentDidMount() {
        this.getArticles();
        console.log('start');
    }

    onSelectionChanged = () => {
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

    changeQueryTypes(type) {
        if (this.state.queryTypes.has(type)) {
            this.state.queryTypes.delete(type)
        } else {
            this.state.queryTypes.add(type)
        }
        this.setState({updateOr: !this.state.updateOr})
    }

    // Summarise

    getSummarise = (task_id, interval = 1000) => {
        fetch(variables.API_URL + `/api/summarise_emb?task_id=${task_id}`,
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
                } else if (res.status == 200) {
                    return res.json()
                } else {
                    throw Error(res.statusText)
                }
            })
        .then((data) => {
          this.setState({
            summarise: data.data,
            message: 'Суммаризация прошла успешно'
          });
        })
        .catch((err) => {
          console.log(err);
          this.setState({summarise: null, message: 'Произошла ошибка при суммаризации'});
        });
    }

    createSummariseQuery() {
        fetch(variables.API_URL + '/api/summarise_emb', {
            method: 'POST',
            headers: {
                'Accept': 'application/json',
                'Content-Type': 'application/json;charset=utf-8',
                'Authorization': `Token ${variables.token}`,
            },
            body: JSON.stringify({
                articles: 'some'
            })
        })
            .then((res) => {
                if (res.status == 200) { return res.json() }
                else { throw Error(res.statusText) }
            })
            .then((result) => {
                var task_id = result.data;
                this.setState({message: 'Отправлено на суммаризацию пожайлуста дождитесь ответа'})
                this.getSummarise(task_id);
            })
            .catch((error) => {
                alert('Ошибка')
            })
    }

    // MarkUp article

    getMarkUp = (task_id, interval = 1000) => {
        fetch(variables.API_URL + `/api/markup?task_id=${task_id}`,
          {
            headers: {
                'Content-Type': 'application/json;charset=utf-8',
                'Authorization': `Token ${variables.token}`,
            },
          }
        )
        .then((res) => {
                if (res.status == 202) {
                    setTimeout(() => {
                      return this.getMarkUp(task_id, interval)
                    }, interval);
                } else if (res.status == 200) {
                    return res.json()
                } else {
                    throw Error(res.statusText)
                }
            })
        .then((data) => {
          try {
            this.setState({
                DetailArticle: data.data,
                message: 'Разметка прошла успешно',
                loading: false,
              });
          } catch {
            console.log('access')
          }
        })
        .catch((err) => {
          console.log(err);
          this.setState({message: 'Произошла ошибка при разметке', loading: false,});
        });
    }

    markUpArticle(DetailArticle) {
        this.setState({loading: true})
        fetch(variables.API_URL + '/api/markup', {
            method: 'POST',
            headers: {
                'Accept': 'application/json',
                'Content-Type': 'application/json;charset=utf-8',
                'Authorization': `Token ${variables.token}`,
            },
            body: JSON.stringify({
                article: DetailArticle
            })
        })
        .then((res) => {
                console.log(res.status)
                if (res.ok) {
                    return res.json()
                } else {
                    throw Error(res.statusText)
                }
            })
        .then((result) => {
          var task_id = result.data;
          this.setState({message: 'Отправлено на суммаризацию пожайлуста дождитесь ответа'})
          this.getMarkUp(task_id);
        })
        .catch((err) => {
          console.log(err);
          this.setState({
            message: 'ошибка при разметке',
            loading: false,
          });
        });
    }

    onRemoveSelected = () => {
        const selectedData = this.gridRef.current.api.getSelectedRows();
        const res = this.gridRef.current.api.applyTransaction({ remove: selectedData });
    }


    render() {
        const {
            token,
            loading,
            query_list,
            articlesInfo,
            articles,
            DetailArticle,
            message,
            summarise,

            queryText,
            queryEndDate,
            queryStartDate,
            queryScore,
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
                      <div className="col-md-2">
                        <div className="row g-0">
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
                                    <h2 class="accordion-header" id="flush-headingThree">
                                      <button class="accordion-button collapsed" data-target='#flush-collapseThree' type="button" data-bs-toggle="collapse" data-bs-target="#flush-collapseThree" aria-expanded="false" aria-controls="flush-collapseThree">
                                        Точность
                                      </button>
                                    </h2>
                                    <div id="flush-collapseThree" className="collapse show multi-collapse" aria-labelledby="flush-headingFour" data-bs-target="#accordionFlushExample">
                                      <div className="accordion-body">
                                        <p>Требуемая точность = {queryScore.toFixed(2)}</p>
                                        <Slider
                                            axis="x"
                                            x={queryScore}
                                            xmax={1}
                                            xmin={0}
                                            xstep={0.01}
                                            onChange={({x}) => this.setState({queryScore: x})}
                                        />
                                      </div>
                                    </div>
                                  </div>
                                  <div className="accordion-item">
                                    <h2 class="accordion-header" id="flush-headingThree">
                                      <button class="accordion-button collapsed" data-target='#flush-collapseThree' type="button" data-bs-toggle="collapse" data-bs-target="#flush-collapseFive" aria-expanded="false" aria-controls="flush-collapseThree">
                                        Выделенные сущности
                                      </button>
                                    </h2>
                                    <div id="flush-collapseFive" className="collapse show multi-collapse" aria-labelledby="flush-headingFour" data-bs-target="#accordionFlushExample">
                                      <div className="accordion-body">
                                        {Object.entries(obj_color).map(tag => <p class="pb-2 mb-3 border-bottom" style={{color: `${tag[1]}`}}>{tag[0]}.</p>)}
                                      </div>
                                    </div>
                                  </div>
                                </div>
                            </aside>

                            <section class="col shadow p-4" style={{backgroundColor: "#fff"}}>
                              <div class="accordion accordion-flush" id="accordion">
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
                                    </div>
                                  </div>
                                <div class="accordion-item">
                                  <h2 class="accordion-header" id="">
                                    <button class="accordion-button collapsed" type="button" data-bs-toggle="collapse" data-bs-target="#flush-collapseSeven" aria-expanded="false" aria-controls="flush-collapseSeven">
                                      Запросы. {message}
                                    </button>
                                  </h2>
                                  <div id="flush-collapseSeven" class="collapse multi-collapse" aria-labelledby="flush-headingSeven" data-bs-target="#accordionFlushExample">
                                    <div class="accordion-body">
                                    {query_list?.map((query, index) =>
                                    query.status === 2?
                                        <p class="pb-2 mb-3 border-bottom" style={{color: 'red'}}>{index + 1} - {query.query}.</p>
                                    :query.status === 1?
                                        <p class="pb-2 mb-3 border-bottom" style={{color: 'green'}}>{index + 1} - {query.query}.</p>
                                    :
                                        <p class="pb-2 mb-3 border-bottom" style={{color: 'black'}}>{index + 1} - {query.query}.</p>
                                    )}
                                    </div>
                                    <div class="accordion-body">
                                        <input className="btn btn-primary" type="submit" value="Очистить" onClick={() => this.clearTask()}/>
                                    </div>
                                  </div>
                                </div>
                              </div>
                              <div>
                              <div class="bd-example">
                                <div class="tab-content" id="myTabContent">
                                  <div class="tab-pane fade active show" id="home" role="tabpanel" aria-labelledby="home-tab">
                                    <div class="container-fluid g-0">
                                        <div className="ag-theme-alpine ag-theme-acmecorp" style={{height: 700}}>
                                            <AgGridReact
                                                ref={this.gridRef}
                                                rowData={articles}
                                                columnDefs={articlesInfo}
                                                pagination={true}
                                                rowSelection={'single'}
                                                onSelectionChanged={this.onSelectionChanged}
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
                                            {summarise?
                                            <>
                                                <p>Summarise</p>
                                                <p>{summarise}</p>
                                            </>
                                            :<input className="btn btn-primary" type="submit" value="Суммаризовать" onClick={() => this.createSummariseQuery()}/>}
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
                                          <a href= { DetailArticle.url } class="card-title link-primary text-decoration-none h5" target="_blank"> { DetailArticle.titl } </a>
                                          <p class="card-text">---------------------------------- </p>
                                          <p class="card-text">Авторы :  { DetailArticle.auth } </p>
                                          <p class="card-text">---------------------------------- </p>
                                          <p class="card-text">Аннотация :  </p>
                                          <p class="card-text" dangerouslySetInnerHTML={{__html: markup_text(DetailArticle.tiab, DetailArticle.annotations)}} />
                                          <p class="card-text">---------------------------------- </p>
                                          <p class="card-text"><small class="text-success">Дата публикации : { DetailArticle.pdat } </small></p>
                                          <p class="card-text"><small class="text-success">Издание : { DetailArticle.jour }</small></p>
                                          <p class="card-text"><small class="text-success">Вид публикации : { DetailArticle.pt }</small></p>
                                          <p class="card-text"><small class="text-success">Страна : { DetailArticle.pl } </small></p>
                                          <p class="card-text"><small class="text-success">{ DetailArticle.mesh } </small></p>
                                          {summarise?
                                            <>
                                                <p>Summarise</p>
                                                <p>{summarise}</p>
                                            </>
                                            :loading?
                                                <p>Loading...</p>
                                            :<input className="btn btn-primary" type="submit" value="Разметить" onClick={() => this.markUpArticle(DetailArticle)}/>
                                          }
                                          <input className="btn btn-primary" type="submit" value="Удалить" onClick={() => this.onRemoveSelected()}/>
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