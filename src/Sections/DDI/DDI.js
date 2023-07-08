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

import { variables } from '../Variables.js';


export class DDIReview extends Component {

    constructor(props) {
        super(props);

        this.gridRef = createRef();
        this.gridAnaliseRef = createRef();
        this.state = {
            token: variables.token,

            //queries
            query_list: [],
            task_id: null,
            articles: [],
            articlesInfo: [
                {field: 'pmid', filter: 'agTextColumnFilter'},
                {field: 'score', filter: 'agTextColumnFilter'},
                {field: 'text', filter: 'agTextColumnFilter'},
            ],

            // Filters
            queryText: 'covid-19',
            queryStartDate: '2022-01-01',
            queryEndDate: null,
            queryTypes: new Set(),
        }
    }

    getArticles = (task_id, interval = 1000) => {
      fetch(variables.API_URL + `/api/ddi_review?task_id=${task_id}`,
            {
                headers: {
                    'Content-Type': 'application/json;charset=utf-8',
                    'Authorization': `Token ${variables.token}`,
                },
            }
          )
          .then(response => {
                console.log(response.status);
                if (response.status == 202) {
                    this.setState({loading: true})
                    setTimeout(() => {
                      return this.getArticles(task_id, interval)
                    }, interval);
                }
                if (response.status == 200) {
                    return response.json()
                } else {
                    throw Error(response.statusText)
                }
          })
          .then(data => {
            this.setState({
                articles: data, DetailArticle: data[0], loading: false
            });
          })
          .catch(error => {
            console.log(error);
                this.setState({ articles: [], DetailArticle: null, loading: false });
          })
    }

    createTask() {
        // Отправляем запрос на сервер для получения статей
        fetch(variables.API_URL + '/api/ddi_review', {
            method: 'POST',
            headers: {
                'Accept': 'application/json',
                'Content-Type': 'application/json;charset=utf-8',
                'Authorization': `Token ${variables.token}`,
            },
            body: JSON.stringify({
                query: this.state.queryText,
            })
            })
            .then(response => response.json())
            .then(data => {
                this.state.query_list.push(this.state.queryText)
                this.setState({
                    task_id: data.data
                });
                alert("You query in queue? please wait to get result");
                this.getArticles(data.data)
            })
            .catch(error => {
                console.log(error);
                this.setState({ task: null });
            }
        )
    }

    componentDidMount() {
        this.getArticles(null);
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

    render() {
        const {
            token,
            query_list,
            articlesInfo,
            articles,
            DetailArticle,

            queryText,
            queryEndDate,
            queryStartDate,
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
                                      Запросы
                                    </button>
                                  </h2>
                                  <div id="flush-collapseSeven" class="collapse multi-collapse" aria-labelledby="flush-headingSeven" data-bs-target="#accordionFlushExample">
                                    <div class="accordion-body">
                                    {query_list?.map(query => <p class="pb-2 mb-3 border-bottom">{query}.</p>)}
                                    </div>
                                  </div>
                                </div>
                              </div>
                              <div>
                              <div class="bd-example">
                                <div class="tab-content" id="myTabContent">
                                  <div class="tab-pane fade active show" id="home" role="tabpanel" aria-labelledby="home-tab">
                                    <div class="container-fluid g-0">
                                        <div className="ag-theme-alpine" style={{height: 700}}>
                                            <AgGridReact
                                                ref={this.gridRef}
                                                rowData={articles}
                                                columnDefs={articlesInfo}
                                                pagination={true}
                                                rowSelection={'single'}
                                                onSelectionChanged={this.onSelectionChanged}
                                            >
                                            </AgGridReact>
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
                                          <a href= { DetailArticle.url } class="card-title link-primary text-decoration-none h5"> { DetailArticle.pmid } </a>
                                          <p class="card-text">---------------------------------- </p>
                                          <p class="card-text">Аннотация :  </p>
                                          <p class="card-text"> { DetailArticle.text } </p>
                                          <p class="card-text">---------------------------------- </p>
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